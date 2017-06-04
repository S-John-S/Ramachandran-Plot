#
# Requires rama500-*.data files from Lovell et al 2003, and an input
# file of phi/psi angles.
#

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

ramachandran.plot <- function(x.scatter, y.scatter,
    x.grid = seq(0, 1, len = nrow(z)), y.grid = seq(0, 1, len = ncol(z)), z.grid,
    xlim = range(x.grid, finite = TRUE), ylim = range(y.grid, finite = TRUE),
    zlim = range(z.grid, finite = TRUE), levels = pretty(zlim, nlevels),
    nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) -
        1), plot.title="", plot.axes, key.title, key.axes, asp = NA,
    xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes,
    ...)
{
    if (missing(z.grid)) {
        stop("no 'z.grid' matrix specified")
    }
    else if (is.list(x.grid)) {
        y.grid <- x.grid$y
        x.grid <- x.grid$x
    }
    if (any(diff(x.grid) <= 0) || any(diff(y.grid) <= 0))
        stop("increasing 'x.grid' and 'y.grid' values expected")

    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)


    if (!is.matrix(z.grid) || nrow(z.grid) <= 1 || ncol(z.grid) <= 1)
        stop("no proper 'z.grid' matrix specified")
    if (!is.double(z.grid))
        storage.mode(z.grid) <- "double"
    filled.contour(as.double(x.grid), as.double(y.grid), z.grid, as.double(levels),
        col = col)

    if (!(missing(x.scatter)) && !(missing(y.scatter))) {
        plot.xy(xy.coords(x.scatter,y.scatter,NULL,NULL,NULL,NULL),
                xlim=xlim, ylim=ylim, xlab="", ylab="", asp=asp,
                type="p", pch=20, cex=0.1)
    }
        
    if (missing(plot.axes)) {
        if (axes) {
            title(main=plot.title, xlab=expression(phi), ylab=expression(psi))
            axis(1, at=c(-180,-90,0,90,180))
            axis(2, at=c(-180,-90,0,90,180))
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
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

#You should edit these filenames to match your own system:

grid.dir <- "/home/john/Downloads/New Folder/top500-angles/pct/rama/"
scatter.filename <- "/home/john/Downloads/New Folder/1HMP.tsv"

scatter.data <- load.scatter(scatter.filename)

png(filename="/home/john/Downloads/New Folder/1HMP.png",width=410, height=410)
#pdf("C:/temp/ramachandran/new/1HMP.pdf")


#Split plot into quadrants:
par(mfrow=c(2,2))

for(rama.type in names(grid.filenames)) {
    col.name = grid.columnnames[rama.type]
    scatter.phi <- scatter.data[which(scatter.data[,"type"]==col.name),"phi"]
    scatter.psi <- scatter.data[which(scatter.data[,"type"]==col.name),"psi"]

    grid.filename <- paste(grid.dir,grid.filenames[rama.type],sep="")
    grid <- load.grid(grid.filename, mid.points)

    #Use small margins to make the plots nice and big, which as a side
    #effect means squeezing the axes labels a bit, and specify a
    #SQUARE plot area (to go with aspect ratio, asp=1)
    par(mar=c(3,3,3,3), mgp=c(1.75,0.75,0), pty="s")

    ramachandran.plot(scatter.phi, scatter.psi,
             x.grid=mid.points, y.grid=mid.points, z.grid=grid,
             plot.title=grid.captions[rama.type],
             levels=grid.levels[rama.type,],
             col=grid.colors[rama.type,])

}

dev.off()
