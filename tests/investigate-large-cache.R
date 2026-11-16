# Optional scale check: writes about 1.3 GiB temporarily.
# Rscript --vanilla tests/investigate-large-cache.R
local({
suppressPackageStartupMessages({library(DBI);library(RSQLite)})
is_assignment <- function(x,name)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name))
definition <- Filter(function(x)is_assignment(x,'server'),as.list(parse('server.R')))[[1]]
statements <- as.list(definition[[3]][[3]])[-1]
tasks <- new.env(parent=globalenv())
for(expression in statements)if(is.call(expression)&&identical(expression[[1]],as.name('<-'))&&
  is.call(expression[[3]])&&identical(expression[[3]][[1]],as.name('function')))eval(expression,tasks)
folder <- tempfile('large-saved-spectra-');dir.create(folder)
on.exit(unlink(folder,recursive=TRUE),add=TRUE)
times <- seq(1.2,82.2,by=1)
scans <- lapply(seq_along(times),function(i)cbind(c(50.0002,1050.0002),c(i,2*i)))
raw <- tasks$investigate_memory_raw(scans,times,'Large.mz5')
cat('Building',length(times),'slices on 1,000,001 bins\n');flush.console()
build <- system.time(cache <- tasks$investigate_cache_build(raw,folder))
cat('Cache files MiB:',file.info(cache$paths)$size/1024^2,'; build seconds:',build[['elapsed']],'\n');flush.console()
for(statistic in c('max','min')) {
  timing <- system.time(result <- tasks$investigate_cache_spectrum(cache,statistic,c(1.2,82.2)))
  scale <- if(statistic=='max')length(times) else 1
  expected <- scale*(1+(result$data$mz-50.0005)/1000)
  stopifnot(max(abs(result$data$intensity-expected))<1e-7)
  cat('PASS:',statistic,'full-run seconds:',timing[['elapsed']],'; blocks:',result$blocks_read,'\n');flush.console()
}
timing <- system.time(region <- tasks$investigate_cache_spectrum(cache,'max',c(1.2,82.2),c(162,165)))
stopifnot(region$blocks_read<result$blocks_read,max(abs(region$data$intensity-length(times)*(1+(region$data$mz-50.0005)/1000)))<1e-7)
cat('PASS: narrow-region seconds:',timing[['elapsed']],'; blocks:',region$blocks_read,'\n');flush.console()
reopened <- tasks$investigate_cache_open(folder,strict=TRUE)
stopifnot(!is.null(reopened))
cat('PASS: completed large caches reopen without mz5 data\n')
if (file.exists('/proc/self/status')) cat(readLines('/proc/self/status')[grepl('^VmHWM:',readLines('/proc/self/status'))], '\n')
unlink(folder,recursive=TRUE)

})
