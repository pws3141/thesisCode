# plotting help: taking from Jonty Rougier blog post

## sort out the graphics windows

graphics.off()

dev.new(width = 5.5, height = 3) # in inches
acrossthetop <- dev.cur()

## plus other sizes ...

## this figure goes across the top

dev.set(acrossthetop)

## draw the figure here

## sort out par

par(mar = c(4, 4, 3, 1), mgp = c(2, 0.7, 0), cex.lab = 0.8,
  cex.axis = 0.8, cex.main = 1, font.main = 1, las = 1)
op <- par(no.readonly = TRUE)

par(op) # reset par values

## then the local par options, then draw the figure

## intercept dev.print

my.dev.print <- function(name, prefix = "foo1") {
  if (PROD) {
    fname <- paste0(prefix, "_", name, ".png")
    dev.print(png, fname,
      width = par("din")[1], height = par("din")[2],
      units = "in", pointsize = 12, res = 144,
      bg = "transparent")
    message("Figure written to ", fname)
  }
  invisible(NULL)
}
