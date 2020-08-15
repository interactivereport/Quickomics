# 3D scatterplots and point identification via rgl

# checked in 23 December 2009 by J. Fox
# 5 January 2010: fixed axis labeling in scatter3d.formula. J. Fox
# 13 May 2010: changed default id.n to conform to showLabels
# 30 July 2010: checks for rgl
# 23 October 2010: added surface.alpha and ellipsoid.alpha arguments
# 2012-03-02: fixed some argument abbreviations. J. Fox
# 2013-02-20: fixed error message, docs for surface.col argument. J. Fox
# 2013-08-20: changed rgl:::rgl.projection to rgl::rgl.projection; more such fixes to come. J. Fox
# 2013-08-31: rgl functions used now exported; got rid of ::: and ::. J. Fox
# 2014-08-04: changed name of identify3d() to Identify3d(). J. Fox
# 2014-08-17: added calls to requireNamespace and :: as needed. J. Fox
# 2014-09-04: J. Fox: empty groups produce warning rather than error
# 2015-12-12: Added axis.ticks argument and code contributed by David Winsemius to add tick labels to axes. J. Fox
# 2016-02-06: Changed call to rgl.clear() to next3d() for compatibility with embedding in HTML. J. Fox
# 2017-06-27: introduced id argument replacing several arguments. J. Fox
# 2017-11-30: use carPalette(), avoid red and green. J. Fox

nice <- function(x, direction=c("round", "down", "up"), lead.digits=1){
	direction <- match.arg(direction)
	if (length(x) > 1) return(sapply(x, nice, direction=direction, lead.digits=lead.digits))
	if (x == 0) return(0)
	power.10 <- floor(log(abs(x),10))
	if (lead.digits > 1) power.10 <- power.10 - lead.digits + 1
	lead.digit <- switch(direction,
		round=round(abs(x)/10^power.10),
		down=floor(abs(x)/10^power.10),
		up=ceiling(abs(x)/10^power.10))
	sign(x)*lead.digit*10^power.10
}


isFALSE <- function(x) length(x) == 1 && is.logical(x) && !isTRUE(x)

applyDefaults <- function(args, defaults, type=""){
    if (isFALSE(args)) return(FALSE)
    names <- names(args)
    names <- names[names != ""]
    if (!isTRUE(args) && !is.null(args) && length(names) != length(args)) warning("unnamed ", type, " arguments, will be ignored")
    if (isTRUE(args) || is.null(names)) defaults
    else defaults[names] <- args[names]
    as.list(defaults)
}


scatter3d <- function(x, ...){
    if (!requireNamespace("rgl")) stop("rgl package missing")
    UseMethod("scatter3d")
}

scatter3d.formula <- function(formula, data, subset, radius, xlab, ylab, zlab, id=FALSE, ...){
    na.save <- options(na.action=na.omit)
    on.exit(options(na.save))
    m <- match.call(expand.dots=FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    m$na.action <- na.pass
    m$id <- m$xlab <- m$ylab <- m$zlab <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    formula <- as.character(formula)
    formula <- paste(formula[2], formula[1], formula[3])
    formula <- as.formula(sub("\\|", "+", formula))
    m$formula <- formula
    X <- eval(m, parent.frame())
    if ("(radius)" %in% names(X)){
        radius <- X[, "(radius)"]
        X <- X[, names(X) != "(radius)"]
    }
    else radius <- 1
    names <- names(X)
    id <- applyDefaults(id, defaults=list(method="mahal", n=2,
        labels=rownames(X), offset = ((100/nrow(X))^(1/3))*0.02), type="id")
    if (missing(xlab)) xlab <- names[2]
    if (missing(ylab)) ylab <- names[1]
    if (missing(zlab)) zlab <- names[3]
    if (ncol(X) == 3) 
        scatter3d(X[,2], X[,1], X[,3], xlab=xlab, ylab=ylab, zlab=zlab, radius=radius, id=id, ...)
    else if (ncol(X) == 4) 
        scatter3d(X[,2], X[,1], X[,3], groups=X[,4], xlab=xlab, ylab=ylab, zlab=zlab, radius=radius, id=id, ...)
    else stop("incorrect scatter3d formula")
}

scatter3d.default <- function(x, y, z,
    xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
    zlab=deparse(substitute(z)), axis.scales=TRUE, axis.ticks=FALSE,
    revolutions=0, bg.col=c("white", "black"),
    axis.col=if (bg.col == "white") c("darkmagenta", "black", "darkcyan")
    else c("darkmagenta", "white", "darkcyan"),
    surface.col=carPalette()[-1],
    surface.alpha=0.5,
    neg.res.col="magenta", pos.res.col="cyan",
    square.col=if (bg.col == "white") "black" else "gray", point.col="yellow",
    text.col=axis.col, grid.col=if (bg.col == "white") "black" else "gray",
    fogtype=c("exp2", "linear", "exp", "none"),
    residuals=(length(fit) == 1), surface=TRUE, fill=TRUE, grid=TRUE, grid.lines=26,
    df.smooth=NULL, df.additive=NULL,
    sphere.size=1, radius=1, threshold=0.01, speed=1, fov=60, 
    fit="linear", groups=NULL, parallel=TRUE, ellipsoid=FALSE, level=0.5, ellipsoid.alpha=0.1,
    # id.method=c("mahal", "xz", "y", "xyz", "identify", "none"), 
    # id.n=if (id.method == "identify") Inf else 0,
    # labels=as.character(seq(along=x)), offset = ((100/length(x))^(1/3)) * 0.02,
    id=FALSE, model.summary=FALSE, ...){
    if (!requireNamespace("rgl")) stop("rgl package missing")
    if (!requireNamespace("mgcv")) stop("mgcv package missing")
  
    id <- applyDefaults(id, defaults=list(method="mahal", n=2,
        labels=as.character(seq(along=x)), offset = ((100/length(x))^(1/3))*0.02), type="id")
    if (isFALSE(id)){
      id.n <- 0
      id.method <- "mahal"
      labels <- NULL
    }
    else{
      labels <- id$labels
      id.method <- id$method
      id.n <- if ("identify" %in% id.method) Inf else id$n
      offset <- id$offset
    }
  
  
#    id.method <- match.arg(id.method)
    if (residuals == "squares"){
        residuals <- TRUE
        squares <- TRUE
    }
    else squares <- FALSE
    summaries <- list()
    if ((!is.null(groups)) && (nlevels(groups) > length(surface.col)))
        stop(sprintf("Number of groups (%d) exceeds number of colors (%d)",
            nlevels(groups), length(surface.col)))
    if ((!is.null(groups)) && (!is.factor(groups))) stop("groups variable must be a factor")
    counts <- table(groups)
    if (any(counts == 0)){
        levels <- levels(groups)
        warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
        groups <- factor(groups, levels=levels[counts != 0])
    }
    bg.col <- match.arg(bg.col)
    fogtype <- match.arg(fogtype)
    if ((length(fit) > 1) && residuals && surface)
        stop("cannot plot both multiple surfaces and residuals")
    xlab  # cause these arguments to be evaluated
    ylab
    zlab
    rgl::next3d()
    rgl::rgl.viewpoint(fov=fov)
    rgl::rgl.bg(color=bg.col, fogtype=fogtype)
    if (id.method == "identify"){
        xg <- x
        yg <- y
        zg <- z
        ggroups <- groups
        glabels <- labels
    }
    valid <- if (is.null(groups)) complete.cases(x, y, z)
    else complete.cases(x, y, z, groups)
    x <- x[valid]
    y <- y[valid]
    z <- z[valid]
    labels <- labels[valid]
    minx <- min(x)
    maxx <- max(x)
    miny <- min(y)
    maxy <- max(y)
    minz <- min(z)
    maxz <- max(z)
    if (axis.scales){
        lab.min.x <- nice(minx)
        lab.max.x <- nice(maxx)
        lab.min.y <- nice(miny)
        lab.max.y <- nice(maxy)
        lab.min.z <- nice(minz)
        lab.max.z <- nice(maxz)
        minx <- min(lab.min.x, minx)
        maxx <- max(lab.max.x, maxx)
        miny <- min(lab.min.y, miny)
        maxy <- max(lab.max.y, maxy)
        minz <- min(lab.min.z, minz)
        maxz <- max(lab.max.z, maxz)
        min.x <- (lab.min.x - minx)/(maxx - minx)
        max.x <- (lab.max.x - minx)/(maxx - minx)
        min.y <- (lab.min.y - miny)/(maxy - miny)
        max.y <- (lab.max.y - miny)/(maxy - miny)
        min.z <- (lab.min.z - minz)/(maxz - minz)
        max.z <- (lab.max.z - minz)/(maxz - minz)
        if (axis.ticks){
            if (axis.scales) {
                x.labels <-  seq(lab.min.x, lab.max.x, 
                    by=diff(range(lab.min.x, lab.max.x))/4)
                x.at <- seq(min.x, max.x, by=nice(diff(range(min.x, max.x))/4))
                rgl::rgl.texts(x.at, -0.05, 0, x.labels, col = axis.col[1])
                
                z.labels <-  seq(lab.min.z, lab.max.z, 
                    by=diff(range(lab.min.z, lab.max.z))/4)
                z.at <- seq(min.z, max.z, by=diff(range(min.z, max.z))/4)
                rgl::rgl.texts(0, -0.1, z.at, z.labels, col = axis.col[3])
                
                y.labels <-  seq(lab.min.y, lab.max.y, 
                    by=diff(range(lab.min.y, lab.max.y))/4)
                y.at <- seq(min.y, max.y, by=diff(range(min.y, max.y))/4)
                rgl::rgl.texts(-0.05, y.at, -0.05, y.labels, col = axis.col[2])
            }
        }
        else {
            rgl::rgl.texts(min.x, -0.05, 0, lab.min.x, col=axis.col[1])
            rgl::rgl.texts(max.x, -0.05, 0, lab.max.x, col=axis.col[1])
            rgl::rgl.texts(0, -0.1, min.z, lab.min.z, col=axis.col[3])
            rgl::rgl.texts(0, -0.1, max.z, lab.max.z, col=axis.col[3])
            rgl::rgl.texts(-0.05, min.y, -0.05, lab.min.y, col=axis.col[2])
            rgl::rgl.texts(-0.05, max.y, -0.05, lab.max.y, col=axis.col[2])
        }
    }
    if (!is.null(groups)) groups <- groups[valid]
    x <- (x - minx)/(maxx - minx)
    y <- (y - miny)/(maxy - miny)
    z <- (z - minz)/(maxz - minz)
    size <- sphere.size*((100/length(x))^(1/3))*0.015
    radius <- radius/median(radius)
    if (is.null(groups)){
        if (size > threshold) rgl::rgl.spheres(x, y, z, color=point.col, radius=size*radius)
        else rgl::rgl.points(x, y, z, color=point.col)
    }
    else {
        if (size > threshold) rgl::rgl.spheres(x, y, z, color=surface.col[as.numeric(groups)], radius=size*radius)
        else rgl::rgl.points(x, y, z, color=surface.col[as.numeric(groups)])
    }
    if (!axis.scales) axis.col[1] <- axis.col[3] <- axis.col[2]
    rgl::rgl.lines(c(0,1), c(0,0), c(0,0), color=axis.col[1])
    rgl::rgl.lines(c(0,0), c(0,1), c(0,0), color=axis.col[2])
    rgl::rgl.lines(c(0,0), c(0,0), c(0,1), color=axis.col[3])
    rgl::rgl.texts(1, 0, 0, xlab, adj=1, color=axis.col[1])
    rgl::rgl.texts(0, 1.05, 0, ylab, adj=1, color=axis.col[2])
    rgl::rgl.texts(0, 0, 1, zlab, adj=1, color=axis.col[3])
    # if (axis.scales){
    #     rgl::rgl.texts(min.x, -0.05, 0, lab.min.x, col=axis.col[1])
    #     rgl::rgl.texts(max.x, -0.05, 0, lab.max.x, col=axis.col[1])
    #     rgl::rgl.texts(0, -0.1, min.z, lab.min.z, col=axis.col[3])
    #     rgl::rgl.texts(0, -0.1, max.z, lab.max.z, col=axis.col[3])
    #     rgl::rgl.texts(-0.05, min.y, -0.05, lab.min.y, col=axis.col[2])
    #     rgl::rgl.texts(-0.05, max.y, -0.05, lab.max.y, col=axis.col[2])
    # }
    if (ellipsoid) {
        dfn <- 3
        if (is.null(groups)){
            dfd <- length(x) - 1
            ell.radius <- sqrt(dfn * qf(level, dfn, dfd))
            ellips <- ellipsoid(center=c(mean(x), mean(y), mean(z)),
                shape=cov(cbind(x,y,z)), radius=ell.radius)
            if (fill) rgl::shade3d(ellips, col=surface.col[1], alpha=ellipsoid.alpha, lit=FALSE)
            if (grid) rgl::wire3d(ellips, col=surface.col[1], lit=FALSE)
        }
        else{
            levs <- levels(groups)
            for (j in 1:length(levs)){
                group <- levs[j]
                select.obs <- groups == group
                xx <- x[select.obs]
                yy <- y[select.obs]
                zz <- z[select.obs]
                dfd <- length(xx) - 1
                ell.radius <- sqrt(dfn * qf(level, dfn, dfd))
                ellips <- ellipsoid(center=c(mean(xx), mean(yy), mean(zz)),
                    shape=cov(cbind(xx,yy,zz)), radius=ell.radius)
                if (fill) rgl::shade3d(ellips, col=surface.col[j], alpha=ellipsoid.alpha, lit=FALSE)
                if (grid) rgl::wire3d(ellips, col=surface.col[j], lit=FALSE)
                coords <- ellips$vb[, which.max(ellips$vb[1,])]
                if (!surface) rgl::rgl.texts(coords[1] + 0.05, coords[2], coords[3], group,
                    col=surface.col[j])
            }
        }
    }
    if (surface){
        vals <- seq(0, 1, length.out=grid.lines)
        dat <- expand.grid(x=vals, z=vals)
        for (i in 1:length(fit)){
            f <- match.arg(fit[i], c("linear", "quadratic", "smooth", "additive"))
            if (is.null(groups)){
                mod <- switch(f,
                    linear = lm(y ~ x + z),
                    quadratic = lm(y ~ (x + z)^2 + I(x^2) + I(z^2)),
                    smooth = if (is.null(df.smooth)) mgcv::gam(y ~ s(x, z))
                    else mgcv::gam(y ~ s(x, z, fx=TRUE, k=df.smooth)),
                    additive = if (is.null(df.additive)) mgcv::gam(y ~ s(x) + s(z))
                    else mgcv::gam(y ~ s(x, fx=TRUE, k=df.additive[1]+1) +
                            s(z, fx=TRUE, k=(rev(df.additive+1)[1]+1)))
                )
                if (model.summary) summaries[[f]] <- summary(mod)
                yhat <- matrix(predict(mod, newdata=dat), grid.lines, grid.lines)
                if (fill) rgl::rgl.surface(vals, vals, yhat, color=surface.col[i], alpha=surface.alpha, lit=FALSE)
                if(grid) rgl::rgl.surface(vals, vals, yhat, color=if (fill) grid.col
                    else surface.col[i], alpha=surface.alpha, lit=FALSE, front="lines", back="lines")
                if (residuals){
                    n <- length(y)
                    fitted <- fitted(mod)
                    colors <- ifelse(residuals(mod) > 0, pos.res.col, neg.res.col)
                    rgl::rgl.lines(as.vector(rbind(x,x)), as.vector(rbind(y,fitted)), as.vector(rbind(z,z)),
                        color=as.vector(rbind(colors,colors)))
                    if (squares){
                        res <- y - fitted
                        xx <- as.vector(rbind(x, x, x + res, x + res))
                        yy <- as.vector(rbind(y, fitted, fitted, y))
                        zz <- as.vector(rbind(z, z, z, z))
                        rgl::rgl.quads(xx, yy, zz, color=square.col, alpha=surface.alpha, lit=FALSE)
                        rgl::rgl.lines(xx, yy, zz, color=square.col)
                    }
                }
            }
            else{
                if (parallel){
                    mod <- switch(f,
                        linear = lm(y ~ x + z + groups),
                        quadratic = lm(y ~ (x + z)^2 + I(x^2) + I(z^2) + groups),
                        smooth = if (is.null(df.smooth)) mgcv::gam(y ~ s(x, z) + groups)
                        else mgcv::gam(y ~ s(x, z, fx=TRUE, k=df.smooth) + groups),
                        additive = if (is.null(df.additive)) mgcv::gam(y ~ s(x) + s(z) + groups)
                        else mgcv::gam(y ~ s(x, fx=TRUE, k=df.additive[1]+1) +
                                s(z, fx=TRUE, k=(rev(df.additive+1)[1]+1)) + groups)
                    )
                    if (model.summary) summaries[[f]] <- summary(mod)
                    levs <- levels(groups)
                    for (j in 1:length(levs)){
                        group <- levs[j]
                        select.obs <- groups == group
                        yhat <- matrix(predict(mod, newdata=cbind(dat, groups=group)), grid.lines, grid.lines)
                        if (fill) rgl::rgl.surface(vals, vals, yhat, color=surface.col[j], alpha=surface.alpha, lit=FALSE)
                        if (grid) rgl::rgl.surface(vals, vals, yhat, color=if (fill) grid.col
                            else surface.col[j], alpha=surface.alpha, lit=FALSE, front="lines", back="lines")
                        rgl::rgl.texts(1, predict(mod, newdata=data.frame(x=1, z=1, groups=group)), 1,
                            paste(group, " "), adj=1, color=surface.col[j])
                        if (residuals){
                            yy <- y[select.obs]
                            xx <- x[select.obs]
                            zz <- z[select.obs]
                            fitted <- fitted(mod)[select.obs]
                            res <- yy - fitted
                            rgl::rgl.lines(as.vector(rbind(xx,xx)), as.vector(rbind(yy,fitted)), as.vector(rbind(zz,zz)),
                                col=surface.col[j])
                            if (squares) {
                                xxx <- as.vector(rbind(xx, xx, xx + res, xx + res))
                                yyy <- as.vector(rbind(yy, fitted, fitted, yy))
                                zzz <- as.vector(rbind(zz, zz, zz, zz))
                                rgl::rgl.quads(xxx, yyy, zzz, color=surface.col[j], alpha=surface.alpha, lit=FALSE)
                                rgl::rgl.lines(xxx, yyy, zzz, color=surface.col[j])
                            }
                        }
                    }
                }
                else {
                    levs <- levels(groups)
                    for (j in 1:length(levs)){
                        group <- levs[j]
                        select.obs <- groups == group
                        mod <- switch(f,
                            linear = lm(y ~ x + z, subset=select.obs),
                            quadratic = lm(y ~ (x + z)^2 + I(x^2) + I(z^2), subset=select.obs),
                            smooth = if (is.null(df.smooth)) mgcv::gam(y ~ s(x, z), subset=select.obs)
                            else mgcv::gam(y ~ s(x, z, fx=TRUE, k=df.smooth), subset=select.obs),
                            additive = if (is.null(df.additive)) mgcv::gam(y ~ s(x) + s(z), subset=select.obs)
                            else mgcv::gam(y ~ s(x, fx=TRUE, k=df.additive[1]+1) +
                                    s(z, fx=TRUE, k=(rev(df.additive+1)[1]+1)), subset=select.obs)
                        )
                        if (model.summary) summaries[[paste(f, ".", group, sep="")]] <- summary(mod)
                        yhat <- matrix(predict(mod, newdata=dat), grid.lines, grid.lines)
                        if (fill) rgl::rgl.surface(vals, vals, yhat, color=surface.col[j], alpha=surface.alpha, lit=FALSE)
                        if (grid) rgl::rgl.surface(vals, vals, yhat, color=if (fill) grid.col
                            else surface.col[j], alpha=surface.alpha, lit=FALSE, front="lines", back="lines")
                        rgl::rgl.texts(1, predict(mod, newdata=data.frame(x=1, z=1, groups=group)), 1,
                            paste(group, " "), adj=1, color=surface.col[j])
                        if (residuals){
                            yy <- y[select.obs]
                            xx <- x[select.obs]
                            zz <- z[select.obs]
                            fitted <- fitted(mod)
                            res <- yy - fitted
                            rgl::rgl.lines(as.vector(rbind(xx,xx)), as.vector(rbind(yy,fitted)), as.vector(rbind(zz,zz)),
                                col=surface.col[j])
                            if (squares) {
                                xxx <- as.vector(rbind(xx, xx, xx + res, xx + res))
                                yyy <- as.vector(rbind(yy, fitted, fitted, yy))
                                zzz <- as.vector(rbind(zz, zz, zz, zz))
                                rgl::rgl.quads(xxx, yyy, zzz, color=surface.col[j], alpha=surface.alpha, lit=FALSE)
                                rgl::rgl.lines(xxx, yyy, zzz, color=surface.col[j])
                            }
                        }
                    }
                }
            }
        }
    }
    else levs <- levels(groups)
    if (id.method == "identify"){
        Identify3d(xg, yg, zg, axis.scales=axis.scales, groups=ggroups, labels=glabels, 
            col=surface.col, offset=offset)
    }
    else if (id.method != "none"){
        if (is.null(groups)) 
            showLabels3d(x, y, z, labels, id.method=id.method, id.n=id.n, col=surface.col[1])
        else {
            for (j in 1:length(levs)){
                group <- levs[j]
                select.obs <- groups == group
                showLabels3d(x[select.obs], y[select.obs], z[select.obs], labels[select.obs], 
                    id.method=id.method, id.n=id.n, col=surface.col[j])
            }
        }
    }
    if (revolutions > 0) {
        for (i in 1:revolutions){
            for (angle in seq(1, 360, length.out=360/speed)) rgl::rgl.viewpoint(-angle, fov=fov)
        }
    }
    if (model.summary) return(summaries) else return(invisible(NULL))
}

# the following function is a slight modification of rgl.select3d() in the rgl package,
#   altered to pass through arguments (via ...) to rgl.select()

car.select3d <- function (...) {
    if (!requireNamespace("rgl")) stop("rgl package is missing")
    rgl::.check3d()
    rect <- rgl::rgl.select(...)
    llx <- rect[1]
    lly <- rect[2]
    urx <- rect[3]
    ury <- rect[4]
    if (llx > urx) {
        temp <- llx
        llx <- urx
        urx <- temp
    }
    if (lly > ury) {
        temp <- lly
        lly <- ury
        ury <- temp
    }
    proj <- rgl::rgl.projection()
    function(x, y, z) {
        pixel <- rgl::rgl.user2window(x, y, z, projection = proj)
        apply(pixel, 1, function(p) (llx <= p[1]) && (p[1] <=
                urx) && (lly <= p[2]) && (p[2] <= ury) && (0 <= p[3]) &&
                (p[3] <= 1))
    }
}

Identify3d  <- function (x, y, z, axis.scales=TRUE, groups = NULL, labels = 1:length(x),
    col = c("blue", "green", "orange", "magenta", "cyan", "red", "yellow", "gray"),
    offset = ((100/length(x))^(1/3)) * 0.02){
    if (!requireNamespace("rgl")) stop("rgl package is missing")
    if (!is.null(groups)){
        counts <- table(groups)
        if (any(counts == 0)){
            levels <- levels(groups)
            warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
            groups <- factor(groups, levels=levels[counts != 0])
        }
    }
    valid <- if (is.null(groups))
        complete.cases(x, y, z)
    else complete.cases(x, y, z, groups)
    labels <- labels[valid]
    x <- x[valid]
    y <- y[valid]
    z <- z[valid]
    groups <- groups[valid]
    minx <- min(x)
    maxx <- max(x)
    miny <- min(y)
    maxy <- max(y)
    minz <- min(z)
    maxz <- max(z)
    if (axis.scales){
        lab.min.x <- nice(minx)
        lab.max.x <- nice(maxx)
        lab.min.y <- nice(miny)
        lab.max.y <- nice(maxy)
        lab.min.z <- nice(minz)
        lab.max.z <- nice(maxz)
        minx <- min(lab.min.x, minx)
        maxx <- max(lab.max.x, maxx)
        miny <- min(lab.min.y, miny)
        maxy <- max(lab.max.y, maxy)
        minz <- min(lab.min.z, minz)
        maxz <- max(lab.max.z, maxz)
        min.x <- (lab.min.x - minx)/(maxx - minx)
        max.x <- (lab.max.x - minx)/(maxx - minx)
        min.y <- (lab.min.y - miny)/(maxy - miny)
        max.y <- (lab.max.y - miny)/(maxy - miny)
        min.z <- (lab.min.z - minz)/(maxz - minz)
        max.z <- (lab.max.z - minz)/(maxz - minz)
    }
    x <- (x - minx)/(maxx - minx)
    y <- (y - miny)/(maxy - miny)
    z <- (z - minz)/(maxz - minz)
    rgl::rgl.bringtotop()
    identified <- character(0)
    groups <- if (!is.null(groups))
        as.numeric(groups[valid])
    else rep(1, length(x))
    repeat {
        f <- car.select3d(button="right")
        which <- f(x, y, z)
        if (!any(which))
            break
        rgl::rgl.texts(x[which], y[which] + offset, z[which], labels[which],
            color = col[groups][which])
        identified <- c(identified, labels[which])
    }
    unique(identified)
}

showLabels3d <- function(x, y, z, labels,
    id.method = "identify", id.n=length(x), col=c("blue"), 
    res=y - mean(y), range.x=range(x), range.z=range(z), 
    offset = ((100/length(x))^(1/3)) * 0.02) {
    if (!requireNamespace("rgl")) stop("rgl package is missing")
    if (id.method == "none") return(NULL)
    if(id.n > 0L) {
        if (missing(labels))
            labels <- as.character(seq(along=x))
        getPoints <- function(w) {
            names(w) <- labels
            iid <- seq(length=id.n)
            ws <- w[order(-w)[iid]]
            match(names(ws), labels)
        }
        ind <-  switch(id.method,
            xz = getPoints(rowSums(qr.Q(qr(cbind(1, x, z))) ^ 2)),
            y = getPoints(abs(res)),
            xyz = union(getPoints(abs(x - mean(x))), union(abs(z - mean(z)),
                getPoints(abs(res)))),
            mahal= getPoints(rowSums(qr.Q(qr(cbind(1, x, y, z))) ^ 2)))
        rgl::rgl.texts(x[ind], y[ind] + offset, z[ind], labels[ind],
            color = col)
        return(labels[ind])
    } 
}

ellipsoid <- function(center=c(0, 0, 0), radius=1, shape=diag(3), n=30){
    if (!requireNamespace("rgl")) "rgl package is missing"
    # adapted from the shapes3d demo in the rgl package
    degvec <- seq(0, 2*pi, length.out=n)
    ecoord2 <- function(p) c(cos(p[1])*sin(p[2]), sin(p[1])*sin(p[2]), cos(p[2]))
    v <- t(apply(expand.grid(degvec,degvec), 1, ecoord2))
    v <- center + radius * t(v %*% chol(shape))
    v <- rbind(v, rep(1,ncol(v)))
    e <- expand.grid(1:(n-1), 1:n)
    i1 <- apply(e, 1, function(z) z[1] + n*(z[2] - 1))
    i2 <- i1 + 1
    i3 <- (i1 + n - 1) %% n^2 + 1
    i4 <- (i2 + n - 1) %% n^2 + 1
    i <- rbind(i1, i2, i4, i3)
    rgl::qmesh3d(v, i)
}
