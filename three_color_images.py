import argparse
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import math
import aplpy

def make_3color_images(field_name,
                       dir="./",
                       ext="drz.fits",
                       red_filter="F814W",
                       green_filter="F475W",
                       blue_filter="F336W",
    ):
    """
    field_name: str
        Name of field.

    dir: str
        Directory to drizzle files.

    ext: str
        Extension for drizzle files. Will assume there are "_" inbetween field
        name, filters, and extension.
    """

    filters = {'red':red_filter, 'green':green_filter, 'blue':blue_filter}

    # check to see that all the required drz files exist
    for i in range(len(filters)):
        drz_file = dir + field_name + "_" + list(filters.values())[i] + "_" + ext
        if not os.path.isfile(drz_file):
            print(drz_file + " was not found")
            break

    # have to rewrite fits files to put science product as the zeroth extension
    # could potentially remove this step if aplpy.make_rgb_cube
    # ever decides to add an hdu input variable
    for i in range(len(filters)):
        drz_file = dir + field_name + "_" + list(filters.values())[i] + "_" + ext
        data = fits.open(drz_file)
        hdu = fits.PrimaryHDU(data["SCI"].data, header=data["SCI"].header)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dir + "temp_" + list(filters)[i] + ".fits", overwrite=True)

    print("Constructing color cube...")
    cubename = dir + field_name + "_colorcube.fits"
    aplpy.make_rgb_cube([dir + "temp_red.fits",
                         dir + "temp_green.fits",
                         dir + "temp_blue.fits"],
                         cubename)

    print("Creating RGB image...")
    imagename = dir + field_name + "_image.png"
    aplpy.make_rgb_image(cubename,
                        imagename,
                        vmin_r=0.,
                        vmin_g=0.,
                        vmin_b=0.,
                        pmax_r=97.,
                        pmax_g=97.5,
                        pmax_b=97.5,
                       )


    fits2 = cubename.replace(".fits", "_2d.fits")

    if "SMC" in field_name:
        # distance to SMC, reference: Scowcroft et al 2016
        dist = 62 # kpc

    if "LMC" in field_name:
        # distance to LMC, reference: ? :)
        dist = 50 # kpc

    # calculate how big 5 pc is in degrees
    five_pc_deg = (5/(dist*1000))*(180/math.pi)

    print("Plotting RGB image...")
    fig = aplpy.FITSFigure(fits2)
    fig.show_rgb(imagename)

    fig.tick_labels.set_font(size='large')
    fig.axis_labels.set_font(size="x-large", weight='medium')
    fig.add_label(0.5, 1.05, text=field_name, relative=True, size='x-large', layer='title')
    fig.add_grid()
    fig.grid.set_xspacing(0.05)
    fig.grid.set_yspacing(0.02)
    fig.grid.set_color("white")
    fig.grid.set_alpha(0.7)
    fig.grid.set_linewidth(1)

    fig.add_scalebar(five_pc_deg)
    fig.scalebar.set_corner("bottom left")
    fig.scalebar.set_color("white")
    fig.scalebar.set_label("5 pc")
    fig.scalebar.set_linewidth(3)
    fig.scalebar.set_font(size='large')

    fig.save(imagename)
    print("Saving " + imagename)

    # clean up the temporary files that were created
    os.remove(dir + "temp_red.fits")
    os.remove(dir + "temp_green.fits")
    os.remove(dir + "temp_blue.fits")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "field",
        type=str,
        help="Field name",
    )
    parser.add_argument(
        "-dir",
        type=str,
        default="./",
        help="Directory",
    )
    parser.add_argument(
        "-ext",
        type=str,
        default="drz.fits",
        help="Drizzle extension",
    )
    parser.add_argument(
        "-r",
        type=str,
        default="F814W",
        help="Red filter",
    )
    parser.add_argument(
        "-g",
        type=str,
        default="F475W",
        help="Green filter",
    )
    parser.add_argument(
        "-b",
        type=str,
        default="F336W",
        help="Blue filter",
    )

    args = parser.parse_args()

    make_3color_images(
        field_name=args.field,
        dir=args.dir,
        ext=args.ext,
        red_filter=args.r,
        green_filter=args.g,
        blue_filter=args.b,
    )
