# Provides helper functions for some of the R scripts included in this repository.

library(ggplot2)
library(extrafont)
library(scales)
library(gplots)

options(error=traceback)

# From http://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual
base_breaks <- function(n=10) {
    function(x) {
		r = log10(range(x, na.rm=TRUE))
		r[1] = r[1] / 2
		r[2] = r[2] * 2
        floor(axisTicks(log10(range(x, na.rm=TRUE)), log=TRUE, n=n))
    }
}

#' Makes a scatterplot
#'
#' Draws a scatterplot.
#' @param data A two-column matrix or data frame, where the first column contains the x values and the second column contains the y values
#' @param draw_best_fit If true, draw a blue best fit line
#' @param draw_diagonal If true, draw a red line representing the diagonal
#' @param draw_diagonal If true, indicate the Pearson correlation between the two variables below the x axis label
#' @keywords graphics, scatterplot
#' @export
scatterplot = function(data, draw_best_fit=TRUE, draw_diagonal=FALSE, show_pearson_correlation=FALSE) {

	if (show_pearson_correlation) {
		pc = paste("(Pearson correlation = ", correlation, ")", sep="")
	} else {
		pc = ""
	}

	plot(x=data[,1], y=data[,2], main="", sub=pc, xlab=colnames(data)[1], ylab=colnames(data)[2])
	fit = lsfit(as.vector(data[,1]), as.vector(data[,2])) # Do a least squares fit.

	if (draw_best_fit) {
		correlation = sprintf("%.3f", cor(data[,1], data[,2]))
		abline(fit, col="blue")	# Color it blue.
	}

	if (draw_diagonal) {
		abline(a=0, b=1, col="red") # Draw a red diagonal.
	}
}

#' Draws a scatterplot using ggplot
#'
#' Draws a scatterplot using ggplot
#' @param data A data frame
#' @param xcol The name of the column for the x axis
#' @param ycol The name of the column for the y axis
#' @param outfile The name of the output file
#' @param bgfill The background color of the plot
#' @param font The font to use in the plot
#' @param pgmx X major gridlines
#' @param pgmy Y major gridlines
#' @keywords graphics, scatterplot
#' @import ggplot2
#' @import extrafont
#' @export
scatterplot_ggplot = function(data, xcol, ycol, outfile, bgfill="white", font="Arial", pgmx = element_blank(), pgmy=element_blank()) {
	xcol = .ggplot_string_backtick(xcol)
	ycol = .ggplot_string_backtick(ycol)
	p = ggplot(data, aes_string(x=xcol, y=ycol)) +
		geom_point() +
		theme(
	        panel.background=element_rect(fill = bgfill), # Make background color white
			panel.border=element_rect(color="black", fill=NA),
	        panel.grid.major = element_line(colour = "black", size=0.5), # Make grid lines black and thin
	      	legend.title=element_blank(), # Remove the legend title
			panel.grid.major.x = pgmx,
			panel.grid.major.y = pgmy,
			text=element_text(family=font)) # Set the font

	.make_image(outfile, plot=p)
}

#' Formats an integer so that it has thousands commas
#'
#' Formats an integer so that it has thousands commas
#' @param x The integer to format
#' @keywords format, number, integer
#' @export
int_with_commas = function (x)
{
	formatC(x, format="d", big.mark=',')
}

#' Draws a horizontal line plot
#'
#' Draws a horizontal line plot
#' @param data Data to plot
#' @param categories_col The name of the column to use as the categories
#' @param values_col The name of the column to use as the values
#' @param outfile The name of the output file
#' @param bgfill The background color of the plot
#' @param font The font to use in the plot
#' @keywords graphics, bar plot
#' @import ggplot2
#' @import extrafont
#' @import gplots
#' @export
horizontal_line_graph = function(data,
	category_col,
	value_col,
	color_col,
	outfile,
	bgfill="white",
	font="Arial",
	x_axis_title=TRUE,
	y_axis_title=FALSE,
    y_axis_title_size=FALSE,
	text_size=10,
	legend_pos="right",
	grid_color="grey90",
	image_height=NULL,
    image_width=NULL,
	scale_color_manual=scale_colour_manual(values = c("#fb6a4a","#FFFF00","#73c475")), #fb6a4a = red; #FFFF00 = yellow; #73c475 = green
	xlogscale=FALSE,
	x_text_angle=90,
	x_axis_vjust=0.5,
	x_text_size=10,
	facet_group="None",
	y_text_size=10,
	vline_intercepts=NULL,
	vlinedf=NULL,
	num_x_ticks=10,
	loglimits=c(NA,NA),
    title=NA
	) {

	px = pretty(c(data[,value_col], vline_intercepts))

	if (xlogscale) {

		scale_x = scale_x_continuous(trans=log_trans(), breaks=base_breaks(n=num_x_ticks), limits=loglimits, labels=function(x) format(x, scientific = FALSE))
	} else {
		if (facet_group == "None") {
			if (!is.null(vline_intercepts)) {
				range_upper = max(c(range(px)[2], max(vline_intercepts)))
			} else {
				range_upper = range(px)[2]
			}

			range_lower = range(px)[1]
			scale_x = scale_x_continuous(breaks=px, limits=c(range_lower, range_upper))
		} else {
			scale_x = scale_x_continuous(breaks=pretty_breaks(), limits=c(0,NA))
		}
	}

	value_col = .ggplot_string_backtick(value_col)
	category_col = .ggplot_string_backtick(category_col)
	color_col = .ggplot_string_backtick(color_col)

	if (x_axis_title) {
		x_axis_title = element_text()
	} else {
		x_axis_title = element_blank()
	}

	if (y_axis_title) {
        if (y_axis_title_size) {
		    y_axis_title = element_text(size=y_axis_title_size)
        } else {
            y_axis_title = element_text()
        }
	} else {
		y_axis_title = element_blank()
	}

	p = ggplot(data, aes_string(x=value_col, y=category_col, color=color_col, group=category_col), position=position_dodge(height=0.1,width=0.1)) +
		scale_x +
		scale_y_discrete() +
	    theme(
			text=element_text(size=text_size, family=font, color="black"), # This sets the text size and font for all elements; can also adjust text size for individual elements if desired
	        panel.background=element_rect(fill = bgfill), # Set background color
	        panel.grid.major=element_line(colour=grid_color, size=0.5), # Make grid lines black and thin
			panel.border=element_rect(color="black", fill=NA),
			legend.position=legend_pos,
			axis.ticks.y=element_blank(), # Remove the y ticks
	        axis.ticks.x=element_blank(), # Remove the x ticks
	    	panel.grid.major.y=element_blank(), # Suppress horizontal grid lines
	    	axis.title.x=x_axis_title, # Set the x axis title,
			axis.title.y=y_axis_title, # Set the x axis title,
	      	legend.title=element_blank(), # Remove the legend title
			legend.key=element_rect(colour=NA, fill=NA),
			axis.text.x=element_text(angle=x_text_angle, size=x_text_size, vjust=x_axis_vjust),
			axis.text.y=element_text(size=y_text_size),
			strip.text.x = element_text(size = 8)
			) +
		geom_line(color="black") + geom_point()

		if (facet_group != "None") {
			p = p + facet_wrap(facet_group, scales="free_x")
			p = p + theme(panel.spacing = unit(1, "lines"))
			if (!is.null(vlinedf)) {
				p = p + geom_vline(aes_string(xintercept=.ggplot_string_backtick("1+ method")), data=vlinedf, color="#FFFF00", alpha=0.7)
                p = p + geom_vline(aes_string(xintercept=.ggplot_string_backtick("2+ methods")), data=vlinedf, color="#73c475", alpha=0.7)
			}
		}

        if (!is.na(title)) {
            p = p + ggtitle(title)
        }

		if (!is.null(scale_color_manual)) {
			p = p + scale_color_manual
		}

		if (!is.null(vline_intercepts)) {
            if (facet_group == "None") {
                num_facets = 1
            } else {
                num_facets = length(unique(data[,facet_group]))
            }

            p = p + geom_vline(xintercept=vline_intercepts, color=rep(c("#FFFF00", "#73c475"), num_facets), alpha=rep(0.7, num_facets*2))

		}

		.make_image(outfile, plot=p, height=image_height, width=image_width)
}




#' Draws a horizontal line plot
#'
#' Draws a horizontal line plot
#' @param data Data to plot
#' @param categories_col The name of the column to use as the categories
#' @param values_col The name of the column to use as the values
#' @param outfile The name of the output file
#' @param bgfill The background color of the plot
#' @param font The font to use in the plot
#' @keywords graphics, bar plot
#' @import ggplot2
#' @import extrafont
#' @import gplots
#' @export
horizontal_line_graph2 = function(data,
	category_col,
	value_col,
	color_col,
	outfile,
	bgfill="white",
	font="Arial",
	x_axis_title=TRUE,
	y_axis_title=FALSE,
    y_axis_title_size=FALSE,
	text_size=10,
	legend_pos="right",
	grid_color="grey90",
	image_height=2,
	scale_color_manual=scale_colour_manual(values = c("#fb6a4a","#73c475")), #fb6a4a = red; #FFFF00 = yellow; #73c475 = green
	xlogscale=FALSE,
	x_text_angle=90,
	x_axis_vjust=0.5,
	x_text_size=10,
	facet_group="None",
	y_text_size=10,
	vline_intercepts=NULL,
	vlinedf=NULL,
	num_x_ticks=10,
	loglimits=c(NA,NA)
	) {

	px = pretty(c(data[,value_col], vline_intercepts))

	if (xlogscale) {

		scale_x = scale_x_continuous(trans=log_trans(), breaks=c(500,1000,2000,5000,8000), limits=c(300,8000))
	} else {
		if (facet_group == "None") {
			if (!is.null(vline_intercepts)) {
				range_upper = max(c(range(px)[2], max(vline_intercepts)))
			} else {
				range_upper = range(px)[2]
			}

			range_lower = range(px)[1]
			scale_x = scale_x_continuous(breaks=px, limits=c(range_lower, range_upper))
		} else {
			scale_x = scale_x_continuous(breaks=pretty_breaks())
		}
	}

	value_col = .ggplot_string_backtick(value_col)
	category_col = .ggplot_string_backtick(category_col)
	color_col = .ggplot_string_backtick(color_col)

	if (x_axis_title) {
		x_axis_title = element_text()
	} else {
		x_axis_title = element_blank()
	}

	if (y_axis_title) {
        if (y_axis_title_size) {
		    y_axis_title = element_text(size=y_axis_title_size)
        } else {
            y_axis_title = element_text()
        }
	} else {
		y_axis_title = element_blank()
	}

	p = ggplot(data, aes_string(x=value_col, y=category_col, color=color_col, group=category_col), position=position_dodge(height=0.1,width=0.1)) +
		scale_x +
		scale_y_discrete() +
	    theme(
			text=element_text(size=text_size, family=font, color="black"), # This sets the text size and font for all elements; can also adjust text size for individual elements if desired
	        panel.background=element_rect(fill = bgfill), # Set background color
	        panel.grid.major=element_line(colour=grid_color, size=0.5), # Make grid lines black and thin
			panel.border=element_rect(color="black", fill=NA),
			legend.position=legend_pos,
			axis.ticks.y=element_blank(), # Remove the y ticks
	        axis.ticks.x=element_blank(), # Remove the x ticks
	    	panel.grid.major.y=element_blank(), # Suppress horizontal grid lines
	    	axis.title.x=x_axis_title, # Set the x axis title,
			axis.title.y=y_axis_title, # Set the x axis title,
	      	legend.title=element_blank(), # Remove the legend title
			legend.key=element_rect(colour=NA, fill=NA),
			axis.text.x=element_text(angle=x_text_angle, size=x_text_size, vjust=x_axis_vjust),
			axis.text.y=element_text(size=y_text_size),
			strip.text.x = element_text(size = 8)
			) +
		geom_line(color="black") + geom_point()

		if (facet_group != "None") {
			p = p + facet_wrap(facet_group, scales="free_x")
			p = p + theme(panel.spacing = unit(1, "lines"))
			if (!is.null(vlinedf)) {
				p = p + geom_vline(aes_string(xintercept=.ggplot_string_backtick("In benchmark")), data=vlinedf, color="#FFFF00", alpha=0.7)
			}
		}

		if (!is.null(scale_color_manual)) {
			p = p + scale_color_manual
		}

		if (!is.null(vline_intercepts)) {
            if (facet_group == "None") {
                num_facets = 1
            } else {
                num_facets = length(unique(data[,facet_group]))
            }

            p = p + geom_vline(xintercept=vline_intercepts, color=rep(c("#73c475"), num_facets), alpha=rep(0.7, num_facets))

		}

		.make_image(outfile, plot=p, height=image_height)
}




#' Draws a bar graph
#'
#' Draws a bar graph
#' @param data A "melted" data frame (using the "melt" function from the reshape2 module)
#' @param xcol The name of the column to use as the bars
#' @param ycol The name of the column to use as the values
#' @param outfile The name of the output file
#' @param bgfill The background color of the plot
#' @param font The font to use in the plot
#' @keywords graphics, bar plot
#' @import ggplot2
#' @import extrafont
#' @import gplots
#' @export
bar_graph = function(data,
	xcol,
	ycol,
	outfile,
	bgfill="white",
	font="Arial",
	fill=NULL,
	colourtype="discrete",
	x_axis_title=TRUE,
	text_size=17,
	legend_pos="right",
	extra=NULL,
	curve=NULL,
	segment=NULL,
	x_axis_hjust=1,
	x_axis_vjust=0,
	position="dodge",
	barlabel=NULL,
	grid_color="grey90",
	x_text_angle=90,
    x_text_size=NA,
	manual_colors=c("#006077", "#0062D0"),
    facet_group=NULL,
    conditional_bold_x_labels=FALSE,
    x_axis_margin1=0,
    x_axis_margin2=0,
    x_axis_margin3=0,
    x_axis_margin4=0,
    continuous_x_breaks=NULL,
    continuous_y_breaks=NULL,
    ylim=c(NA,NA),
    xlim=c(NA,NA),
    height=NULL,
    width=NULL,
    x_grid_lines=FALSE
) {

    if (is.na(x_text_size)) {
        x_text_size = text_size
    }

	xcol = .ggplot_string_backtick(xcol)
	ycol = .ggplot_string_backtick(ycol)

    if (x_grid_lines) {
        x_grid_lines = element_line()
    } else {
        x_grid_lines = element_blank()
    }

	if (x_axis_title) {
		x_axis_title = element_text()
	} else {
		x_axis_title = element_blank()
	}

	if (is.null(fill)) {
		gb = geom_bar(fill="black", position=position, stat="identity")
	} else {
		fill = .ggplot_string_backtick(fill)
		gb = geom_bar(color="black", mapping=aes_string(fill=fill), position=position, stat="identity")
	}

	p = ggplot(data, aes_string(x=xcol, y=ycol)) +
	    gb +
	    theme(
			text=element_text(size=text_size, family=font), # This sets the text size and font for all elements; can also adjust text size for individual elements if desired
	        panel.background=element_rect(fill = bgfill), # Set background color
	        panel.grid.major=element_line(colour=grid_color, size=0.5), # Make grid lines grey and thin
			legend.position=legend_pos,
			axis.ticks.y=element_blank(), # Remove the y ticks
	        #axis.ticks.x=element_blank(), # Remove the x ticks
	    	panel.grid.major.x=x_grid_lines, # Suppress vertical grid lines
	    	axis.title.x=x_axis_title, # Set the x axis title
	      	legend.title=element_blank(), # Remove the legend conditional_bold_x_labels
			axis.text.x=element_text(angle=x_text_angle, hjust=x_axis_hjust, vjust=x_axis_vjust, size=x_text_size, margin=margin(x_axis_margin1,x_axis_margin2,x_axis_margin3,x_axis_margin4)),
            panel.border=element_rect(color="black", fill=NA),
            strip.background = element_blank(),
            )

        if (conditional_bold_x_labels) {
            p = p + theme(axis.text.x=element_text(face=ifelse(data$bold_x=="Yes","bold.italic","plain")))
        }

		if (colourtype == "continuous") {
			p = p + scale_fill_gradientn(colours=colorpanel(30, "green", "red4", "red"))
		} else {
			p = p + scale_fill_manual(values=manual_colors)
		}

		if (!is.null(facet_group)) {
			p = p + facet_wrap(facet_group)
		}

		if (!is.null(extra)) {
			p = p + extra
		}

		if (!is.null(curve)) {
			p = p + curve
		}

		if (!is.null(segment)) {
			p = p + segment
		}

		if (!is.null(barlabel)) {
			p = p + geom_text(aes_string(label=barlabel))
		}

        if (!is.null(continuous_x_breaks)) {
            p = p + scale_x_continuous(breaks=continuous_x_breaks, limits=xlim)
        }

        if (!is.null(continuous_y_breaks)) {
            p = p + scale_y_continuous(breaks=continuous_y_breaks, limits=ylim)
        }

		.make_image(outfile, plot=p, height=height, width=width)
}

#' Draws a line plot using ggplot
#' @keywords graphics, line plot
#' @import scales
#' @import ggplot2
#' @import extrafont
#' @export
density_plot = function(
    data,
    column,
    fill_col,
    outfile,
    height=NULL,
    width=NULL,
    alpha=0.1,
    ylab=NULL,
    xmin=NULL,
    xmax=NULL
    ) {
        column = .ggplot_string_backtick(column)
        fill_col = .ggplot_string_backtick(fill_col)
        p = ggplot(data, aes_string(column, fill=fill_col)) + geom_density(alpha=alpha)

        if (!is.null(ylab)) {
            p = p + ylab(ylab)
        }

        if (!is.null(xmin) || !is.null(xmax)) {
            p = p + xlim(xmin, xmax)
        }

        .make_image(outfile, plot=p, height=height, width=width)
}

#' Draws a line plot using ggplot
#'
#' Draws a line plot using ggplot
#' @param data A data frame
#' @param xcol The name of the column for the x axis
#' @param ycol The name of the column for the y axis
#' @param outfile The name of the output file
#' @param groupby Which variable to group by?
#' @param nxbreaks Number of x breaks
#' @param nybreaks Number of y breaks
#' @param bgfill The background color of the plot
#' @param font The font to use in the plot
#' @param pgmx X major gridlines
#' @param pgmy Y major gridlines
#' @keywords graphics, line plot
#' @import scales
#' @import ggplot2
#' @import extrafont
#' @export
line_or_scatter_plot = function(
	data,
	xcol,
	ycol,
	outfile,
	groupby=NULL,
    point_shapes=NULL,
    xbreaks=NULL,
	nxbreaks=10,
	nybreaks=10,
	bgfill="white",
	font="Arial",
	text_size=15,
    ylab=NULL,
	pgmx=element_line(colour = "lightgrey", size=0.15),
	pgmy=element_line(colour = "lightgrey", size=0.15),
	box_x_min=1,
	box_x_max=0,
	rectfillcolor="dodgerblue2",
	geom_line=FALSE,
	geom_point=TRUE,
	geom_point_alpha=1, # Lower values make the points somewhat transparent, which is good for when the points overlap
	ylogbreaks=NULL,
	xrev=FALSE,
	ybreaks=NULL,
	manual_colors=NULL,
    vline_intercepts=NULL,
    hline_intercepts=NULL,
    ylim1=NA,
    ylim2=NA,
    xlim1=NA,
    xlim2=NA,
    geom_hexbin=FALSE,
    facet_group="None",
    title=NULL,
    legend_title=FALSE,
    point_size=1,
    geom_label=NULL,
    geom_label_size=1.9,
    height=NULL,
    width=NULL,
    legend_position="right",
    num_facet_cols=NULL,
    label_wrap_gen_width=NULL,
    facet_text_size=15,
    jitter=FALSE
	) {

	if (!is.null(ylogbreaks)) {
		scale_y = scale_y_log10(breaks = ylogbreaks)
	} else {
		scale_y = scale_y_continuous(breaks = pretty_breaks(n=nybreaks), limits=c(ylim1, ylim2))
	}

	if (xrev) {
		scale_x = scale_x_reverse(breaks = pretty_breaks(n=nxbreaks))
	} else {
        if (!is.null(xbreaks)) {
            scale_x = scale_x_continuous(breaks = xbreaks, limits=c(xlim1, xlim2))
        } else {
            scale_x = scale_x_continuous(breaks = pretty_breaks(n=nxbreaks), limits=c(xlim1, xlim2))
        }
	}

	xcol = .ggplot_string_backtick(xcol)
	ycol = .ggplot_string_backtick(ycol)

    if (legend_title) {
        legend_title = element_text()
    } else {
        legend_title = element_blank()
    }

    if (!is.null(groupby)) {
        groupby = .ggplot_string_backtick(groupby)
    }
	#ycolp1 = paste(ycol, "+3", sep="") # Syntax that works for the ribbon. To be potentially implemented later.
	#ycolm1 = paste(ycol, "-3", sep="") # Syntax that works for the ribbon. To be potentially implemented later.

	p = ggplot(data, aes_string(x=xcol, y=ycol, group=groupby)) +
		#geom_ribbon(aes_string(color=groupby, ymin=ycolm1, ymax=ycolp1, fill=groupby), alpha=0.4) + # Ribbon - left out here - to be potentially implemented later
		scale_x +
		scale_y +
		theme(
            panel.spacing = unit(0.7, "lines"),
			text=element_text(size=text_size, family=font),
			panel.background=element_rect(fill = bgfill), # Make background color white
			panel.border=element_rect(color="black", fill=NA),
	      	legend.title=legend_title,
            legend.position=legend_position,
			panel.grid.major.x=pgmx, # Set the y gridlines
			panel.grid.major.y=pgmy, # Set the y gridlines
            plot.margin=margin(10,10,10,10),
			legend.key=element_rect(colour = NA, fill = NA) # Get rid of the gray background on the legend key
			) # Set the font

    if (!is.null(ylab)) {
        p = p + ylab(ylab)
    }

    if (!is.null(title)) {
        p = p + ggtitle(title)
    }

	if (geom_line) {
		p = p + geom_line(aes_string(color=groupby))
	}

    if (geom_hexbin) {
        p = p + geom_hex() + scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(10))
    }

    if (facet_group != "None") {
        p = p + facet_wrap(facet_group, scales="free", ncol=num_facet_cols, labeller=label_wrap_gen(width=label_wrap_gen_width)) + theme(strip.text.x=element_text(size=facet_text_size))
    }

	if (geom_point) {
		p = p + geom_point(aes_string(color=groupby, shape=point_shapes), alpha=geom_point_alpha, size=point_size)

		if (!is.null(manual_colors)) {
			p = p + scale_color_manual(values=manual_colors)
		}

	}

    if (!is.null(vline_intercepts)) {
        p = p + geom_vline(xintercept=vline_intercepts, linetype="dotted")
    }

    if (!is.null(hline_intercepts)) {
        p = p + geom_hline(yintercept=hline_intercepts, linetype="dotted")
    }

    if (!is.null(geom_label)) {
        p = p + geom_label(data=geom_label, aes(x=x, y=y, label=lab), size=geom_label_size, inherit.aes=FALSE, hjust=0)
    }

    if (jitter) {
        p = p + geom_jitter()
    }

# Colour choices for box fill: blue4, deepskyblue
	if (box_x_min < box_x_max) {
		rect = data.frame(xmin=box_x_min, xmax=box_x_max, ymin=-Inf, ymax=Inf)
		p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=rectfillcolor, alpha=0.2, inherit.aes=FALSE)
	}

	.make_image(outfile, plot=p, height=height, width=width)
}


#' @param data A data frame
#' @param xcol The name of the column for the x axis
#' @param ycol The name of the column for the y axis
#' @param outfile The name of the output file
#' @param groupby Which variable to group by?
#' @param nxbreaks Number of x breaks
#' @param nybreaks Number of y breaks
#' @param bgfill The background color of the plot
#' @param font The font to use in the plot
#' @param pgmx X major gridlines
#' @param pgmy Y major gridlines
#' @keywords graphics, line plot
#' @import scales
#' @import ggplot2
#' @import extrafont
#' @export
histogram = function(
	data,
    column,
    outfile=NULL,
    num_bins=10,
    facet_group="None",
    title="",
    font="Arial",
    text_size=15,
    bgfill="white",
    pgmx=element_line(colour = "lightgrey", size=0.15),
    pgmy=element_line(colour = "lightgrey", size=0.15),
    xlab=NA,
    ylab=NA,
    ybreaks=NULL,
    xbreaks=NULL,
    color="black",
    percentage=TRUE,
    x_axis_angle = 0,
    x_axis_hjust = 0.5,
    x_axis_vjust = 0,
    center=NA,
    facet_name_size=8,
    x_axis_text_size=8
	) {

    column = .ggplot_string_backtick(column)

    if (percentage) {
        #percentage_thing = aes(y=..count../sum(..count..))
        percentage_thing = aes(y=100*(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])
    } else {
        percentage_thing = NULL
    }

	p = ggplot(data, aes_string(column)) +
        geom_histogram(percentage_thing, bins=num_bins, fill=color, center=center) +
        ggtitle(title) +
        theme(
    		text=element_text(size=text_size, family=font),
    		panel.background=element_rect(fill = bgfill), # Make background color white
    		panel.border=element_rect(color="black", fill=NA),
    		panel.grid.major.x=pgmx, # Set the y gridlines
    		panel.grid.major.y=pgmy, # Set the y gridlines
            axis.text.x = element_text(angle=x_axis_angle, hjust=x_axis_hjust, vjust=x_axis_vjust, size=x_axis_text_size)
    	)

    if (!is.na(xlab)) {
        p = p + xlab(xlab)
    }

    if (!is.na(ylab)) {
        p = p + ylab(ylab)
    }

    if (!is.null(ybreaks)) {
        p = p + scale_y_continuous(breaks=ybreaks)
    }

    if (!is.null(xbreaks)) {
        p = p + scale_x_continuous(breaks=xbreaks)
    }

    if (facet_group != "None") {
        p = p + facet_wrap(facet_group, scales="free")
        p = p + theme(strip.text.x = element_text(size=facet_name_size))
    }

	.make_image(outfile, plot=p)
}

################################################################################
# PRIVATE FUNCTIONS
################################################################################

# This function is needed because strings with spaces need to be quoted with backticks to be acceptable to ggplot
.ggplot_string_backtick = function(str) {
		return(paste("`", str, "`", sep=""))
}

.make_image = function(outfile, plot, height=NULL, width=NULL) {
	if (!is.null(height) && !is.null(width)) {
		suppressMessages(ggsave(outfile, plot=plot, height=height, width=width))
	} else if (!is.null(height)) {
		suppressMessages(ggsave(outfile, plot=plot, height=height))
	} else if (!is.null(width)) {
		suppressMessages(ggsave(outfile, plot=plot, width=width))
	} else {
		suppressMessages(ggsave(outfile, plot=plot))
	}

	embed_fonts(outfile)
}

.default_breaks = function(a, b) {
		seq(a, b, (a - b)/5)
}
