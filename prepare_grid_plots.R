f = function(t) {((t - 2.1)^3 + (t-5)/3)^2 + 3}
fw = function(t) {(((t - 3.1))^3 + t^2/2.2 - 5*t/2.2)^2 + 3}
fww = function(t) {(t-6)^2 * (t-5)^2 + t^2}
xx = seq(2, 8, length.out=201)



fww1 = function(t) {2*(t-6)*(t-5)^2 + 2*(t-6)^2*(t-5) + 2*t}
fww2 = function(t) {2*(t-5)^2 + 4*(t-6)*(t-5) + 4*(t-6)*(t-5) + 2*(t-6)^2 + 2}

fwwf = function(t, start=1) {fww2(start)/2*(t- (start - fww1(start)/fww2(start)))^2 + fww(start) - fww2(start)/2*(start- (start - fww1(start)/fww2(start)))^2 }


xlim = c(2, 8)
ylim = c(0, 150)

xcoord = function(x) {(x - xlim[[1]]) / diff(xlim)}
ycoord = function(y) {(y - ylim[[1]]) / diff(ylim)}

grid.newpage()
pushViewport(viewport())
drawIt <- function(gp) {
  pushViewport(viewport())
  grid.draw(gp)
  upViewport()
}

gplot <- gTree(
  x = NULL,
  y = NULL,
  childrenvp = vpTree(
    plotViewport(c(5, 4, 4, 2), name = "plotRegion"),
    vpList(
      viewport(
        name="dataRegion",
        xscale = xlim,
        yscale = ylim))),
  children = gList(
    xaxisGrob(vp = "plotRegion::dataRegion"),
    yaxisGrob(vp = "plotRegion::dataRegion")))

gplot <- addGrob(
  gplot,
  linesGrob(x=xcoord(xx),
            y=ycoord(fww(xx)),
            vp = "plotRegion::dataRegion",
            gp = gpar(lty="dotted")))

iterate_gd = function(gp, gamma, x_last) {
  finished = FALSE
  grad = fww1(x_last)
  
  while (!finished) {
    
    gamma = gamma / 0.8
    trial = fww(x_last - grad * gamma)
    if (trial > fww(x_last) + sign(grad)* 0.2*gamma*grad^2) {
      
      gamma = gamma * 0.8^2
    } else {
      x_next = x_last - grad * gamma
      # path = c(path, x_next)
      finished = TRUE
    }
  }
  
  gp <- addGrob(
    gp,
    linesGrob(
      x = xcoord(c(x_last, x_next)),
      y = ycoord(c(fww(x_last), fww(x_last) - gamma*grad^2)),
      gp = gpar(
        col = "purple",
        lty="dotted"),
      vp = "plotRegion::dataRegion"))
  gp = addGrob(
    gp,
    pointsGrob(
      x = x_last,
      y = fww(fww(x_last)),
      gp = gpar(
        col="red"
      ),
      vp = "plotRegion::dataRegion"))
  
  list("gamma" = gamma, "x_new" = x_next, "plot" = gp)
}

# first step
iterate_nr = function(gp, xloc) {
  
  # lines(xx, fwwf(xx, start=xloc), lty=2)
  # xnew = xloc - fww1(xloc)/fww2(xloc)
  # points(xloc, fww(xloc), col='red', pch=16)
  
  gp <- addGrob(
    gp,
    linesGrob(
      x = xcoord(xx),
      y = ycoord(fwwf(xx, start=xloc)),
      gp = gpar(
        col = "purple",
        lty="dotted"),
      vp = "plotRegion::dataRegion"))
  gp = addGrob(
    gp,
    pointsGrob(
      x = xloc,
      y = fww(xloc),
      gp = gpar(
        col="red"
      ),
      vp = "plotRegion::dataRegion"))
  
  
  list(x_new = xloc - fww1(xloc)/fww2(xloc), plot = gp)
}

path = c(2.5)
gamma = 0.005
