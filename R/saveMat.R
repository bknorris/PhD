library(rmatio)
names(mss) <- paste0("layer", seq_along(mss))
write.mat(mss, filename = 'RENAMEstats.mat', compression = FALSE, version = c("MAT5"))