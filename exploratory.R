# Visualización de la información - ITBA
# VAST 2018 - Mini challange 1
# Author: Julián Ailán

rm(list = ls()); gc()

#### Dataset cleanup. ####
inputfile  <- '~/Documents/infoviz/tp/option1/MC1 2018/AllBirdsv4.csv'
outputfile <- '~/Documents/infoviz/tp/option1/own/allbirdsfiltered.csv'
df <- read.csv(inputfile)

# Move every vocalization type to lower case.
df[sapply('\\s*\\wall\\s*$', grepl, df$Vocalization_type), ]$Vocalization_type <- 'call'
df[sapply('Song\\s*$', grepl, df$Vocalization_type), ]$Vocalization_type <- 'song'
df[sapply('\\wall, song\\s*$', grepl, df$Vocalization_type), ]$Vocalization_type <- 'call, song'
df <- df[!(df$Vocalization_type %in% c('scold', 'bill-snapping', 'drumming', '?')), ]

write.csv(df, file = outputfile)
rm(inputfile, outputfile, df)

#### Sound analysis. ####
library(tuneR)
library(signal)
library(seewave)

path  <- '~/Documents/infoviz/tp/option1/MC1 2018/ALL BIRDS'
files <- list.files(path = path, pattern = '*.mp3', full.names = T, recursive = F)
pattern <- '/ALL BIRDS/(.*?)(-)(\\d*).mp3'
points.file <- '~/Documents/infoviz/tp/option1/own/points.csv'

durations <- numeric(length(files))
for(i in 1:length(files)) {
  data <- readMP3(files[i])
  s <- data@left
  fs <- data@samp.rate
  durations[i] <- (length(s) - 1) / fs
}
plot(density(durations), xlab = 'Duration (sec)', ylab = 'Frequency', main = 'Distribution of audio length')

mp3ToPoint <- function(f, file, writetofile, pattern) {
  data <- readMP3(f)
  snd  <- data@left
  fs   <- data@samp.rate
  
  # Remove the DC component/offset of the signal.
  snd <- snd - mean(snd)
  # Normalize to [-1,1] amplitude.
  snd <- snd / (2 ^ (data@bit - 1))
  # Duration in miliseconds
  duration <- 1000 * (0:(length(snd) - 1)) / fs
  # Number of points to use for the FFT. Make sure to use a power of 2.
  n  <- floor(sqrt(length(snd))) ^ 2
  # Apply a high pass filter to remove low frequency noise.
  if(((length(snd) - 1) / fs) < 200) {
    snd <- ffilter(snd, fs, from = 0, to = 500, bandpass = F) 
  }
  # Compute Fast Fourier Transform.
  ft <- fft(snd)
  uniquepts <- ceiling((n + 1) / 2)
  ft <- ft[1:uniquepts]
  freq <- (0:(uniquepts - 1)) * (fs / n)
  # Keep only the real part, normalize so as not to depend on the number of points, and
  # compute the signal power.
  ft <- (abs(ft) / n) ^ 2
  # Move the power to dBs.
  ftdb <- 10 * log10(ft)
  # Point consisting on the frequency (in kHz) where the signal had more power.
  species <- regmatches(f, regexec(pattern, f))[[1]][2]
  recordingid <- regmatches(f, regexec(pattern, f))[[1]][4]
  instance <- data.frame('species' = species, 
                         'power' = max(ftdb), 
                         'frequency' = freq[which(ftdb == max(ftdb))] / 1000,
                         'id' = recordingid)
  if(instance$frequency > 0.5) {
    if(writetofile){
      write.table(instance, file = file, append = T, row.names = F, col.names = F, sep = ',') 
    } 
  }
  return(instance)
}

lapply(files[1:length(files)], mp3ToPoint, points.file, T, pattern)
rm(path, files, pattern, durations, data, fs, i, s, files)

#### Analysis of species by Power and Frequency. ####
inputfile  <- points.file
df <- read.csv(inputfile, header = F)
colnames(df) <- c('species', 'power', 'frequency', 'id')

representation <- data.frame('species' = character(0),
                             'power' = numeric(0),
                             'frequency' = numeric(0),
                             stringsAsFactors = F)
for(species in df$species) {
  p <- mean(df[df$species == species, ]$power)
  f <- mean(df[df$species == species, ]$frequency)
  representation[nrow(representation) + 1, ] <- c(species, p, f)
}
representation <- unique(representation)
representation$power <- as.double(representation$power)
representation$frequency <- as.double(representation$frequency)
rm(inputfile, df, species, f, p)

#### Classification of Kasios' test birds. ####
path    <- '~/Documents/infoviz/tp/option1/MC1 2018/Test Birds from Kasios'
files   <- list.files(path = path, pattern = '*.mp3', full.names = T, recursive = F)
pattern <- '/Test Birds from Kasios/(.*?)(-)(\\d*).mp3'

test.birds <- data.frame('power' = numeric(0),
                         'frequency' = numeric(0))
for(testbird in files) {
  testbirdattr <- mp3ToPoint(testbird, file = points.file, writetofile = T, pattern = pattern)
  test.birds[nrow(test.birds) + 1, ] <- c(testbirdattr$power, testbirdattr$frequency)
}

similar.average.bird <- data.frame('distance' = numeric(0),
                                   'testsp' = character(0),
                                   'repsp' = character(0),
                                   'testfreq' = numeric(0),
                                   'repfreq' = numeric(0),
                                   stringsAsFactors = F)
for(i in 1:nrow(test.birds)) {
  distance <- sqrt(((test.birds$power[i] - representation$power) ^ 2) + 
                   ((test.birds$frequency[i] - representation$frequency) ^ 2))
  min.distance <- min(distance)
  position <- which(distance == min.distance)
  similar.average.bird[nrow(similar.average.bird) + 1, ] <- c(round(min.distance, 5), 
                                                              i,
                                                              representation[position, ]$species,
                                                              round(test.birds$frequency[i], 5),
                                                              round(representation[position, ]$frequency, 5))
}
rm(path, files, pattern, testbirdattr, distance, min.distance, position, i, points.file, testbird)
write.csv(similar.average.bird, 
          file = '~/Documents/infoviz/tp/option1/own/similarities.csv', 
          row.names = F)
