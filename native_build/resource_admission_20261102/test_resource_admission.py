"""Resource-race regressions with fake native execution; no GPU calculations."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import resource_admission as admission

OOM = ('PBEQ CUDA failed: 3 Insufficient GPU memory: need 2.478 GiB '
       'plus 2 GiB reserve, free 3.479 GiB\nERROR STOP 84\n')


class ResourceAdmissionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.workspace = Path(self.temporary.name)
        self.case = self.workspace / 'energy'
        self.case.mkdir()
        (self.workspace / 'status.json').write_text(json.dumps({'stage': 'mutant_pb_energy', 'pid': 123}))
        (self.case / 'pb.inp').write_text('unchanged cold solver deck')
        (self.workspace / 'validated_cache').write_text('previous complete energy')
        self.arguments = ('executable', self.workspace, 'mutant', self.case, 'pb', 1)
        self.addCleanup(patch.stopall)
        patch.object(admission.energy, '_available', return_value=32 * admission.GIB).start()
        patch.object(admission, 'gpu_available', return_value=6 * admission.GIB).start()
        self.sleep = patch.object(admission.time, 'sleep').start()

    def invoke(self, function):
        with patch.object(admission.energy, '_run_process', function):
            admission.install()
            return admission.energy._run_process(*self.arguments)

    def test_detects_exact_gpu_memory_admission_error(self):
        (self.case / 'pb.log').write_text(OOM)
        needed = admission.memory_rejection(RuntimeError('pb: CHARMM exited with code 84'), self.case, 'pb')
        self.assertGreater(needed, (2.478 + 2) * admission.GIB)
        self.assertLess(needed, admission.MIN_GPU)

    def test_code_84_alone_is_not_retryable(self):
        (self.case / 'pb.log').write_text('Solver failed to converge\nERROR STOP 84\n')
        self.assertIsNone(admission.memory_rejection(RuntimeError('pb: CHARMM exited with code 84'), self.case, 'pb'))

    def test_other_exit_code_is_not_retryable(self):
        (self.case / 'pb.log').write_text(OOM)
        self.assertIsNone(admission.memory_rejection(RuntimeError('pb: CHARMM exited with code 83'), self.case, 'pb'))

    def test_incomplete_log_is_not_retryable(self):
        (self.case / 'pb.log').write_text(OOM.replace('ERROR STOP 84', ''))
        self.assertIsNone(admission.memory_rejection(RuntimeError('pb: CHARMM exited with code 84'), self.case, 'pb'))

    def test_vdw_exit_84_is_not_misclassified(self):
        (self.case / 'pb.log').write_text(OOM)
        self.assertIsNone(admission.memory_rejection(RuntimeError('vdw: CHARMM exited with code 84'), self.case, 'vdw'))

    def test_waits_before_launch_and_restores_stage(self):
        calls = []
        stages = []
        def native(*args):
            calls.append(args)
            return {'returncode': 0}
        def observe_sleep(_):
            stages.append(json.loads((self.workspace / 'status.json').read_text())['stage'])
            self.assertEqual(calls, [])
        self.sleep.side_effect = observe_sleep
        with patch.object(admission.energy, '_available', side_effect=[29 * admission.GIB, 32 * admission.GIB, 32 * admission.GIB]), \
             patch.object(admission, 'gpu_available', side_effect=[4 * admission.GIB, 6 * admission.GIB]):
            result = self.invoke(native)
        self.assertEqual(result['returncode'], 0)
        self.assertEqual(len(calls), 1)
        self.assertEqual(stages, ['waiting_native_resources'] * 2)
        state = json.loads((self.workspace / 'status.json').read_text())
        self.assertEqual(state['stage'], 'mutant_pb_energy')
        self.assertEqual(state['pid'], 123)

    def test_oom_retries_same_native_deck_and_preserves_log(self):
        calls = []
        def native(*args):
            calls.append(args)
            if len(calls) == 1:
                (self.case / 'pb.log').write_text(OOM)
                raise RuntimeError('pb: CHARMM exited with code 84')
            (self.case / 'pb.log').write_text('complete, to be audited by original caller')
            return {'returncode': 0}
        self.assertEqual(self.invoke(native), {'returncode': 0})
        self.assertEqual(calls, [self.arguments, self.arguments])
        self.assertEqual((self.case / 'pb.inp').read_text(), 'unchanged cold solver deck')
        self.assertEqual((self.workspace / 'validated_cache').read_text(), 'previous complete energy')
        logs = list(self.case.glob('pb.gpu_oom.*.log'))
        self.assertEqual(len(logs), 1)
        self.assertEqual(logs[0].read_text(), OOM)
        retry = json.loads((self.workspace / 'native_resource_retries.jsonl').read_text())
        self.assertEqual(retry['retry'], 1)
        self.assertEqual(retry['stem'], 'mutant')

    def test_numerical_failure_propagates_without_retry(self):
        calls = []
        def native(*args):
            calls.append(args)
            (self.case / 'pb.log').write_text('Nonfinite energy\nERROR STOP 84\n')
            raise RuntimeError('pb: CHARMM exited with code 84')
        with self.assertRaisesRegex(RuntimeError, 'code 84'):
            self.invoke(native)
        self.assertEqual(len(calls), 1)
        self.assertEqual(list(self.case.glob('pb.gpu_oom.*.log')), [])

    def test_host_admission_race_retries_prelaunch_rejection(self):
        calls = []
        def native(*args):
            calls.append(args)
            if len(calls) == 1:
                raise RuntimeError(admission.HOST_REJECTION)
            return {'returncode': 0}
        self.assertEqual(self.invoke(native), {'returncode': 0})
        self.assertEqual(len(calls), 2)
        self.assertFalse((self.workspace / 'native_resource_retries.jsonl').exists())

    def test_persistent_oom_has_bounded_retries(self):
        calls = []
        def native(*args):
            calls.append(args)
            (self.case / 'pb.log').write_text(OOM)
            raise RuntimeError('pb: CHARMM exited with code 84')
        with self.assertRaisesRegex(RuntimeError, 'code 84'):
            self.invoke(native)
        self.assertEqual(len(calls), admission.MAX_GPU_RETRIES + 1)
        self.assertEqual(len(list(self.case.glob('pb.gpu_oom.*.log'))), admission.MAX_GPU_RETRIES)

    def test_gpu_monitor_failure_blocks_native_launch(self):
        def native(*args):
            self.fail('Native execution must not start with unknown GPU resources')
        with patch.object(admission, 'gpu_available', side_effect=RuntimeError('GPU monitor unavailable')):
            with self.assertRaisesRegex(RuntimeError, 'GPU monitor unavailable'):
                self.invoke(native)

    def test_wait_can_be_interrupted(self):
        def native(*args):
            self.fail('Native execution started before admission')
        self.sleep.side_effect = InterruptedError('interrupted')
        with patch.object(admission.energy, '_available', return_value=29 * admission.GIB):
            with self.assertRaises(InterruptedError):
                self.invoke(native)

    def test_install_is_idempotent(self):
        def native(*args):
            return {'returncode': 0}
        with patch.object(admission.energy, '_run_process', native):
            admission.install()
            once = admission.energy._run_process
            admission.install()
            self.assertIs(admission.energy._run_process, once)
            self.assertEqual(once(*self.arguments), {'returncode': 0})


if __name__ == '__main__':
    unittest.main()
