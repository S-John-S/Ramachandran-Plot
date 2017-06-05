load.scatter <- function(filename) {
  data <- read.table(filename, header=FALSE, comment.char="", row.names=1,
                     colClasses=c("character","numeric","numeric","factor"),
                     col.names=c("name","phi","psi","type"))
  return(data)    
}

load.grid <- function(filename, mid.points) {
  data <- as.matrix(read.table(filename, header=FALSE, comment.char="#",
                               colClasses=c("numeric","numeric","numeric"),
                               col.names=c("phi","psi","value")))
  
  stopifnot(mid.points == seq(-179,179,2))
  stopifnot(min(data[,"phi"])==-179)
  stopifnot(min(data[,"psi"])==-179)
  stopifnot(max(data[,"phi"])==+179)
  stopifnot(max(data[,"psi"])==+179)
  
  # Want to turn data (list of pairs) into grid (matrix)
  grid <- array(0.0, dim = c(length(mid.points),length(mid.points)))
  
  #This works, but is about three times slower than the next solution
  #system.time(for (row.number in 1:nrow(data)) {
  #    i = match(mid.points, data[row.number,"phi"])
  #    j = match(mid.points, data[row.number,"psi"])
  #    grid[i,j] <- data[row.number,"value"]
  #})
  
  #This is fast, but assumes the file is in a precise order.
  row.number <- 0
  for (i in 1:180) for (j in 1:180) {
    row.number <- row.number+1
    if (row.number %% 100 == 0) {
      #Only check every few hundred rows,
      #doing 32,400 checks is too slow!
      stopifnot(data[row.number,"phi"]==mid.points[i])
      stopifnot(data[row.number,"psi"]==mid.points[j])
    }
    grid[i,j] <- data[row.number,"value"]
  }
  remove(row.number)
  remove(data)
  return(grid)
}

mid.points <- seq(-179,179,2)

grid.filenames <- c(General="rama500-general.data",
                    Glycine="rama500-gly-sym.data",
                    Proline="rama500-pro.data",
                    Pre.Pro="rama500-prepro.data")

grid.captions <- c(General="General",
                   Glycine="Glycine (Symmetric)",
                   Proline="Proline",
                   Pre.Pro="Pre-Proline")

grid.columnnames <- c(General="General",
                      Glycine="Glycine",
                      Proline="Proline",
                      Pre.Pro="Pre-Pro")

grid.levels <- t(cbind(General=c(0, 0.0005, 0.02, 1),
                       Glycine=c(0, 0.002,  0.02, 1),
                       Proline=c(0, 0.002,  0.02, 1),
                       Pre.Pro=c(0, 0.002,  0.02, 1)))

grid.colors <- t(cbind(General=c('#FFFFFF','#B3E8FF','#7FD9FF'),
                       Glycine=c('#FFFFFF','#FFE8C5','#FFCC7F'),
                       Proline=c('#FFFFFF','#D0FFC5','#7FFF8C'),
                       Pre.Pro=c('#FFFFFF','#B3E8FF','#7FD9FF')))

grid.dir <- "/home/john/Downloads/New Folder/top500-angles/pct/rama/"

scatter.filename <- "/home/john/Downloads/New Folder/1HMP.tsv"

scatter.data <- load.scatter(scatter.filename)

par(mfrow=c(2,2))

for(rama.type in names(grid.filenames)) {
  col.name = grid.columnnames[rama.type]
  scatter.phi <- scatter.data[which(scatter.data[,"type"]==col.name),"phi"]
  scatter.psi <- scatter.data[which(scatter.data[,"type"]==col.name),"psi"]
  grid.filename <- paste(grid.dir,grid.filenames[rama.type],sep="")
  grid <- load.grid(grid.filename, mid.points)
  par(mar=c(3,3,3,3), mgp=c(1.75,0.75,0), pty="s")
  plot(x=scatter.phi, y=scatter.psi, xlim=c(-180,180), ylim=c(-180,180), main=grid.captions[rama.type],  xlab=expression(phi), ylab=expression(psi), pch=20, cex=0.1, asp=1.0)
  .filled.contour(x=mid.points,y=mid.points,z=grid,levels=grid.levels[rama.type,],col=grid.colors[rama.type,])
  
}


