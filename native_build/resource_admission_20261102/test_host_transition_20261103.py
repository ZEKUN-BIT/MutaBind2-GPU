"""Regression for the Nov 3 VDW-to-PB host-memory admission failure."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import resource_admission as admission


class HostTransitionTests(unittest.TestCase):
    def test_memory_drop_between_vdw_and_pb_waits_and_keeps_vdw(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            case = workspace / 'energy'
            case.mkdir()
            state_path = workspace / 'status.json'
            state_path.write_text(json.dumps({'stage': 'mutant_pb_energy', 'pid': 123}))
            (case / 'vdw.inp').write_text('frozen VDW input')
            (case / 'pb.inp').write_text('frozen PB input')
            calls = []
            observed = []

            def native(executable, root, stem, path, mode, threads):
                calls.append(mode)
                (path / (mode + '.log')).write_text('completed ' + mode)
                return {'returncode': 0}

            def sleeping(_):
                observed.append(json.loads(state_path.read_text())['stage'])
                self.assertEqual(calls, ['vdw'])
                self.assertEqual((case / 'vdw.log').read_text(), 'completed vdw')
                self.assertFalse((case / 'pb.log').exists())

            with patch.object(admission.energy, '_run_process', native), \
                 patch.object(admission.energy, '_available', side_effect=[
                     32 * admission.GIB, 29 * admission.GIB, 28 * admission.GIB, 32 * admission.GIB]), \
                 patch.object(admission, 'gpu_available', return_value=6 * admission.GIB), \
                 patch.object(admission.time, 'sleep', side_effect=sleeping):
                admission.install()
                for mode in ('vdw', 'pb'):
                    admission.energy._run_process('binary', workspace, 'mutant', case, mode, 1)

            self.assertEqual(calls, ['vdw', 'pb'])
            self.assertEqual(observed, ['waiting_native_resources'] * 2)
            self.assertEqual(json.loads(state_path.read_text())['stage'], 'mutant_pb_energy')
            self.assertEqual((case / 'vdw.inp').read_text(), 'frozen VDW input')
            self.assertEqual((case / 'pb.inp').read_text(), 'frozen PB input')
            self.assertEqual((case / 'pb.log').read_text(), 'completed pb')

    def test_unrecognized_host_errors_still_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            case = workspace / 'energy'
            case.mkdir()
            calls = []

            def native(*args):
                calls.append(args)
                raise RuntimeError('Native process exceeded memory limit')

            with patch.object(admission.energy, '_run_process', native), \
                 patch.object(admission.energy, '_available', return_value=32 * admission.GIB), \
                 patch.object(admission, 'gpu_available', return_value=6 * admission.GIB):
                admission.install()
                with self.assertRaisesRegex(RuntimeError, 'exceeded memory limit'):
                    admission.energy._run_process('binary', workspace, 'mutant', case, 'pb', 1)
            self.assertEqual(len(calls), 1)


if __name__ == '__main__':
    unittest.main()
