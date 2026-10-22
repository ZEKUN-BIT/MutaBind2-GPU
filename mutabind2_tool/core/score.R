# Load the original trained R model; receive all filenames as process arguments.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Expected mode, model, feature table, output and contribution flag")
suppressPackageStartupMessages(library(randomForest))
mode <- args[1]
if (!(mode %in% c("single", "multiple"))) stop("Unknown model mode")
model_env <- new.env()
load(args[2], envir = model_env)
model <- model_env[[if (mode == "single") "mutabinds.rf" else "mutabindm.rf"]]
test <- read.delim(args[3], check.names = FALSE, stringsAsFactors = FALSE)
if (!nrow(test)) stop("Feature table contains no rows")
features <- all.vars(delete.response(model$terms))
missing <- setdiff(features, names(test))
if (length(missing)) stop(paste("Missing feature columns:", paste(missing, collapse = ", ")))
for (name in features) {
    value <- suppressWarnings(as.numeric(test[[name]]))
    if (any(!is.finite(value))) stop(paste("Feature contains missing or non-finite values:", name))
    test[[name]] <- value
}
test$SunDDG <- as.numeric(predict(model, test[, features, drop = FALSE]))
if (any(!is.finite(test$SunDDG))) stop("The model returned a non-finite prediction")
if ("interface" %in% names(test)) {
    test$Interface <- test$interface
    test$Confidence <- ifelse(grepl("no", test$interface, fixed = TRUE), "Low", "High")
}
test$Deleterious <- ifelse(abs(test$SunDDG) >= 1.5, "yes", "no")
if (args[5] == "yes") {
    options(rgl.useNULL = TRUE)
    suppressPackageStartupMessages(library(forestFloor))
    set.seed(1)
    train <- model_env$data[, features, drop = FALSE]
    # forestFloor requires matching column classes, including integer counts.
    for (name in features) {
        if (is.integer(train[[name]])) {
            if (any(abs(test[[name]] - round(test[[name]])) > 1e-8)) {
                stop(paste("Expected integer-valued feature:", name))
            }
            test[[name]] <- as.integer(test[[name]])
        }
    }
    ff <- forestFloor(model, train, test[, features, drop = FALSE],
                      bootstrapFC = TRUE, calc_np = TRUE)
    fc <- as.data.frame(tail(ff$FCmatrix, nrow(test)))
    names(fc) <- features
    labels <- c(ddVdw = "DDE_vdw", ddPb = "DDG_solv", Foldx_fold = "DDG_fold",
                FoldX_fold = "DDG_fold", provean = "CS", ACC_wt = "SA_com_wt",
                ACCmp_wt = "SA_part_wt", len_mp_contact_op_wt10 = "N_cont_wt", dVdwWt = "dE_vdw_wt")
    for (name in features) test[[labels[[name]]]] <- round(fc[[name]], 4)
}
write.table(test, args[4], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
