fun.len.arrs <- function(xvel, yvel, prob1 = 25, prob2 = 75, scale.factor = 2) {
        # xvel: wind speed x velocity at arrow points
        # yvel: wind speed y velocity at arrow points
        # prob1: lower quantile prob where arrows start fading to grey
        # prob2: upper quantile prob where arrows start increasing thickness
        # scale.factor: ratio between length of shortest and longest arrows

	stopifnot(length(xvel) == length(yvel))

    xy.vel <- cbind(xvel, yvel)
	len <- nrow(xy.vel)

	scale_by <- apply(xy.vel, 1, function(x) sqrt(sum(x^2)))

    # find lower and upper quantiles
	lower_quant <- quantile(x = scale_by, prob = prob1/100, names = FALSE)
	upper_quant <- quantile(x = scale_by, prob = prob2/100, names = FALSE)

    # find index permutation to put scale_by into increasing order
	vel_norm.order <- order(scale_by)

    # sort scale_by into increasing order
	vel_norm.sort <- scale_by[vel_norm.order]

    # logical vector for finding lower and upper quantiles
	which.lower_quant <- vel_norm.sort < lower_quant
	which.upper_quant <- vel_norm.sort > upper_quant

    # how many arrows in each quantile
	num.lower_quant <- sum(which.lower_quant)
	num.upper_quant <- sum(which.upper_quant)

    # create function that increases arrow size
    # start at relative length '1',
    # gradient
    linear.function <- function(x, x1, x2, y1) {
            a <- (y1 - 1) / (x2 - x1)
            # intercept
            b <- 1 - a * x1
            a * x + b
    }

    # vector of length of each arrow
	len.arrs <- rep(1, length = len)
	col.vec <- rep("gray0", length = len)
	lwd.vec <- rep(1, length = len)
	for(i in 1:len) {
        # obtain position of (size of) arrow
        lower_tmp <- which.lower_quant[i]
        upper_tmp <- which.upper_quant[i]
        # first, consider the arrows in lower quantile
		if (isTRUE(lower_tmp)) {
			len.arrs[i] <- 1 #vel_norm.sort[i]
            # colours get lighter as move towards 0% quantile
			col.vec[i] <- paste("gray", 
                                floor((- 70 / num.lower_quant) * (i - 1) + 70), sep = "")
        } 
        # now, consider the arrows in upper quantile
        if (isTRUE(upper_tmp)) {
			len.arrs[i] <- scale.factor #* vel_norm.sort[i]
            # arrows get thicker as move towards 100% quantile
            x1_tmp <- len - num.upper_quant
            x2_tmp <- len
            y1_tmp <- 2
			lwd.vec[i] <- linear.function(x = i, x1 = x1_tmp, x2 = x2_tmp, y1 = y1_tmp)
        } 
        # now, remaining vectors inbetween the two quantiles
        if (isFALSE(lower_tmp) && isFALSE(upper_tmp)) {
            # start at relative length '1',
            # move to relative length 'scale.factor',
            # over [prob1, prob2] quantile distance
            # move to relative length 'scale.factor',
            # over [prob1, prob2] quantile distance
            x1_tmp <- num.lower_quant
            x2_tmp <- len - num.upper_quant
            y1_tmp <- scale.factor
            len.arrs[i] <- linear.function(x = i, x1 = x1_tmp, x2 = x2_tmp,
                                           y1 = y1_tmp) #* vel_norm.sort[i]
        }
	}
    # order(vel_norm.order) reinstates original ordering
	list(len.arrs[order(vel_norm.order)], col.vec[order(vel_norm.order)], 
         lwd.vec[order(vel_norm.order)])
}

ggplot.arrows <- function(windu, windv, narrx = 15, narry = 15, scale.factor = 2) {
        # windu: windu matrix (i.e. one matrix from windu array)
        # windv: windv matrix (i.e. one matrix from windu array)
        # narrx: number of arrows in x direction
        # narry: number of arrows in y direction
        # arrlen: length of arrows (what units?)
        # lon is MSLP longitude -> cell mid points
        # lat is MSLP latitude -> cell mid points

	u2 <- windu
	v2 <- windv

	#n <- dim(input)[1] / 3
	nx <- dim(u2)[1]
	ny <- dim(v2)[2]

    # split [0, 1] into narrx parts, with half cells at left and right
    arrows.x <- seq(from = 0.5/narrx, to = (narrx - 0.5)/narrx, by = 1/narrx)
    # rep s.t. can pair with arrows.y
	arrows.x <- rep(arrows.x, each = narry)
    # round to determine which cell number the arrows will be in
	nx.arrows.x <- round(nx * arrows.x, digits = 0)

    # split [0, 1] into narry parts, with half cells at bottom and top
    arrows.y <- seq(from = 0.5/narry, to = (narry - 0.5)/narry, by = 1/narry)
    # rep s.t. can pair with arrows.x
	arrows.y <- rep(arrows.y, times = narrx)
    # round to determine which cell number the arrows will be in
	ny.arrows.y <- round(ny * arrows.y, digits = 0)

    # pair the two arrow positions 
    index <- cbind(nx.arrows.x, ny.arrows.y)

    # Find the respective x and y wind speed at the arrow positions
	windu_selected <- u2[index]
	windv_selected <- v2[index]

    winduv_selected <- cbind(windu_selected, windv_selected)
    norm.winduv_selected <- apply(winduv_selected, 1, norm, type = "2")
    # scale s.t. max = 1
    norm.winduv_selected <- norm.winduv_selected / max(norm.winduv_selected)

    # cell number -> longitude and latitude positions
    lon_span <- 3.75
    lon <- seq(from = 1.875, by = 3.75, length = 96)
    for(i in 1:length(lon)) {
    	if(lon[i] > 180) {
    		lon[i] <- lon[i] - 360
    	}
    }
    #lon <- lon - lon_span / 2
    lat_span <- 2.5
    lat <- seq(from = 88.75, by = -lat_span, length = 72)
    #
	lon.x <- lon[nx.arrows.x]
	lat.y <- lat[ny.arrows.y]

	# set modulus of arrows equal to size determined by arrlen
	# scale_by <- (1 / (arrlen/100)) * sqrt(windu_selected^2 + windv_selected^2)
	# differing len or arrs
	arrow_plot.parameters <- fun.len.arrs(xvel = windu_selected, yvel = windv_selected,
                                                    scale.factor = scale.factor)
	scale_by <- arrow_plot.parameters[[1]]
	col.vec <- arrow_plot.parameters[[2]]
	lwd.vec <- arrow_plot.parameters[[3]]

	# rescale arrows
    winduv_selected.scale <- winduv_selected / max(abs(winduv_selected))
	#arrow_length.windu <- winduv_selected.scale * (1 + scale_by * norm.winduv_selected)
	#arrow_length.windv <- windv_selected.scale * (1 + scale_by * norm.winduv_selected)

    # want to keep winduv direction, but have arrow length between [1, scale.factor]
    # set norm of each winduv row to be 1
    winduv_selected.scale <- t(apply(winduv_selected, 1, function(x) x / sqrt(sum(x^2))))
    # scale norm by scale_by
    # '*' multiplies row i by element i
    arrow_xy <- winduv_selected.scale * scale_by

    #scale_vector <- (1 + scale_by * norm.winduv_selected)
    #arrow_length <- winduv_selected * scale_vector
    # scale arrow length s.t.
    #arrow_length.scale

	# set end points
	# note - and +, relates to way `image' plots matrix
	arrow.x1 <- lon.x + arrow_xy[, 1]
	arrow.y1 <- lat.y + arrow_xy[, 2]

	#output
	output <- list(lon.x, lat.y, arrow.x1, arrow.y1, col.vec, lwd.vec)
	names(output) <- c("x0", "y0", "x1", "y1", "col", "lwd")
	output
	# plot image and arrows
	#arrows(x0 = arrows.x, y0 = (1-arrows.y), 
	#	x1 = arrow.x1, y1 = (1-arrow.y1), length = 0.05,
	#	col = col.vec, lwd = lwd.vec)
}
