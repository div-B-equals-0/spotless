import sys
import os.path
import argparse
import numpy as np
import scipy.ndimage as ndi

try:
    # Prefer the astropy version of fits library
    import astropy.io.fits as pyfits
except ImportError:
    # But fall back to pyfits version if not found
    import pyfits

import skimage.filter
import skimage.morphology
from skimage.morphology import erosion, dilation
import pyregion

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Remove cosmic rays and other bad pixels from an image""")

parser.add_argument("fitsfile", type=str, help="""Name of image FITS 
                    file (sans extension)""")
parser.add_argument("--hdu-index", type=int, default=0, 
                    help="""Which HDU to use from the FITS file""")
parser.add_argument("--output-id", type=str, default="", help="""Extra
                    string to add to output filenames to aid in layer
                    identification""")
parser.add_argument("--method", choices=("thresh", "edge", "segment"), default="edge",
                    help="Algorithm to use to find the bad pixels")
parser.add_argument("--onlybadpix", action="store_true", help="""Only
                    calculate the bad pixel map - do not replace
                    pixels in the image""")
parser.add_argument("--threshold", type=float, default=None, 
                    help="Assume any pixel above this level is bad")
parser.add_argument("--dmax", type=int, default=5, 
                    help="""Maximum diameter of features to zap.
                    Leave alone any roughly circular objects that are
                    larger than this. 
                    """)
parser.add_argument("--ddmax", type=int, default=10, 
                    help="""Absolute maximum diameter of features to
                    zap.  Leave alone any objects that are larger than
                    this, whatever their shape may be.
                    """)
parser.add_argument("--data-range", type=float, default=None, 
                    nargs=2, metavar=("MIN", "MAX"),  
                    help="Range for data scaling")
parser.add_argument("--use-log-scale", action="store_true",
                    help="Use logarithmic data scaling")
parser.add_argument("--segment-pars", type=float, default=(0.0, 1.0),
                    nargs=2, metavar=("LO", "HI"), help="""For
                    'segment' method only: thresholds of scaled data
                    to seed the good/bad regions """)
parser.add_argument("--edge-pars", type=float, default=(1.0, 0.1, 0.2),
                    nargs=3, metavar=("SIGMA", "LOW_THRESHOLD",
                                      "HIGH_THRESHOLD"),
                    help="""For 'edge' method only: parameters for the
                    Canny edge detection algorithm.  See:
                    http://scikits-image.org/docs/dev/auto_examples/plot_canny.html
                    """)
parser.add_argument("--thick-edges", action="store_true", help="""Make
                    the edges be 3 pixels wide instaed of the default 1""")
parser.add_argument("--reject-filaments", action="store_true",
                    help="""Try to reject objects that look filamentary, since they are
                    probably not cosmic rays""")
parser.add_argument("--allow-shadows", action="store_true",
                    help="""Also remove objects that are darker than
                    their surroudings (not generally advised,
                    especially if you have dark globules in your
                    image!)""")
parser.add_argument("--clip-negative", action="store_true",
                    help="""Also remove all negative pixels""")
parser.add_argument("--exclude-regions-from-file", type=str, default=None,
                    help="""Read DS9 regions from a file, which are to
                    be marked as definite good pixels""" )
parser.add_argument("--include-regions-from-file", type=str, default=None,
                    help="""Read DS9 regions from a file, which are to
                    be marked as definite bad pixels""" )
parser.add_argument("--verbose", "-v", action="store_true", 
                    help="Print informative progress messages")
parser.add_argument("--debug", "-d", action="store_true", 
                    help="Save auxiliary images of intermediate steps")
parser.add_argument("--multi-hdu", "-m", action="store_true", 
                    help="""Work in multi-HDU mode.  
This assumes that the image is in the "SCI" HDU in the input file (the argument --hdu-index is ignored).  
All additional HDUs in the input file are copied through to the output file.
Only one output file is written, all auxilliary arrays 
("edges", "labels", "badpix", etc) are written as additional HDUs in the same file.
""")


cmd_args = parser.parse_args()

if cmd_args.multi_hdu:
    hdulist = pyfits.open(cmd_args.fitsfile + ".fits")
    inhdu = hdulist["SCI"]
    outprefix = "{fitsfile}-{output_id}".format(**vars(cmd_args))
else:
    inhdu = pyfits.open(cmd_args.fitsfile + ".fits")[cmd_args.hdu_index]
    outprefix = "{fitsfile}-cleaned-{method}{output_id}".format(**vars(cmd_args))

inhdu.header["CRREJECT"] = (True, "Cosmic ray rejection by spotless.py")
inhdu.header["CRMETHOD"] = (cmd_args.method, "CR rejection method")

def save_aux_image(image, label):
    """
    Save an auxiliary image, such as the bad pixel mask.

    If in multi-HDU mode, just append it to the list of HDUs to be written later.
    Otherwise, actually write it to a file.
    """
    if cmd_args.multi_hdu:
        hdulist.append(pyfits.ImageHDU(image, name=label))
    else:
        pyfits.PrimaryHDU(image).writeto(outprefix + "-" + label + ".fits", clobber=True)



if cmd_args.verbose:
    print "Finding bad pixels by the '{method}' method".format(**vars(cmd_args))

## First find the bad pixel mask
if cmd_args.method == "thresh":
    if cmd_args.threshold is None:
        raise ValueError, "The --threshold argument is obligatory with --method thresh"
    badpix = inhdu.data > cmd_args.threshold
else:
    # Remaining methods use scikits image
    # which requires float arrays in range 0.0 -> 1.0
    if not cmd_args.data_range[-1]:
        # fallback to using the actual data limits if no range specified
        datamin, datamax = inhdu.data.min(), inhdu.data.max()
    else:
        datamin, datamax = cmd_args.data_range
    # Scale so that [datamin, datamax] -> [0, 1]
    if cmd_args.use_log_scale:
        # logarithmic scale 
        scaled_data = np.log10(inhdu.data/datamin)/np.log10(datamax/datamin)
    else:
        # linear scale
        scaled_data = (inhdu.data - datamin) / (datamax - datamin)
    # Clip all values outside of the scaled range [0, 1]
    scaled_data[scaled_data < 0.0] = 0.0
    scaled_data[scaled_data > 1.0] = 1.0
    # Any NaNs or Infs also get clamped to the minimum value
    scaled_data[~np.isfinite(scaled_data)] = 0.0
    if cmd_args.verbose:
        print "Data scaled to range [{:.2e}, {:.2e}]".format(datamin, datamax)
    if cmd_args.debug:
        # Save the elevation map to a FITS file for debugging purposes
        save_aux_image(scaled_data, "scaled")

    # 8-neighbor structuring element
    square3 = skimage.morphology.square(3)
    # 4-neighbor structuring element
    cross3 = skimage.morphology.diamond(1)
    if cmd_args.method == "edge":
        # Edge detection followed by filling in of holes
        sig, lo, hi = cmd_args.edge_pars
        # use Canny edge detector
        edges = skimage.filter.canny(scaled_data, sigma=sig, 
                                     low_threshold=lo, high_threshold=hi) 
        if cmd_args.verbose:
            print "Edge detection with Canny method complete"
        if cmd_args.thick_edges:
            # Make the edges 3-pix wide instead of 1
            edges = skimage.morphology.dilation(edges.astype(np.uint8), 
                                                skimage.morphology.square(3)
                                            ).astype(bool)
            if cmd_args.verbose:
                print "Dilation of edges to 3 pixels complete"

        # Save the edges to a FITS file for debugging purposes
        if cmd_args.debug:
            save_aux_image(edges.astype(int), "edges")

        badpix = ndi.binary_fill_holes(edges)
        if cmd_args.verbose:
            print "Filling of holes complete"
        
        # Now remove all lone 1-pixel edges by an erosion followed by a dilation
        badpix = dilation(erosion(badpix.astype(np.uint8), cross3), square3).astype(bool)

    elif cmd_args.method == "segment":
        # Segmentation via watershed transform of elevation_map.  This
        # works by growing regions outwards from seed pixels, using a
        # map of the image gradient magnitude to define basins. 
        #
        # This method does not work very well for cosmic rays  
        #
        # First, we use Sobel filter to find the magnitude of the
        # gradient of the image
        elevation_map = skimage.filter.sobel(scaled_data)
        if cmd_args.verbose:
            print "Elevation map calculation via Sobel gradient filter complete"

        # Save the elevation map to a FITS file for debugging purposes
        if cmd_args.debug:
            save_aux_image(elevation_map, "sobel")

        # Now we mark pixels that are definitely good or bad
        markers = np.zeros_like(scaled_data, dtype=np.uint8)
        loval, hival = cmd_args.segment_pars
        markers[scaled_data < loval] = 1
        markers[scaled_data > hival] = 2
        # Then we segment the image by growing out from the marked
        # values
        segmentation = skimage.morphology.watershed(elevation_map, markers)

        # Save the segmentation to a FITS file for debugging purposes
        if cmd_args.debug:
            save_aux_image(segmentation, "segments")

        if cmd_args.verbose:
            print "Segmentation via watershed filter complete"
        # The rest is the same as with the edge detection method
        badpix = ndi.binary_fill_holes(segmentation-1)
        if cmd_args.verbose:
            print "Filling of holes complete"

    if not cmd_args.threshold is None:
        # Also add into bad pixels all above the threshold
        badpix = badpix | (inhdu.data > cmd_args.threshold)
        if cmd_args.verbose:
            print "All values higher than {:.2e} also added to bad pixels".format(cmd_args.threshold)

    if cmd_args.clip_negative:
        # Also add into bad pixels all negative pixels
        badpix = badpix | (inhdu.data < 0.0)
        if cmd_args.verbose:
            print "All negative values also added to bad pixels"


# reset status from bad -> good in all ring-fenced regions
if cmd_args.exclude_regions_from_file:
    regions = pyregion.open(cmd_args.exclude_regions_from_file)
    mask = regions.get_mask(hdu=inhdu)
    badpix = badpix & ~mask
    if cmd_args.verbose:
        print "Regions from {} removed from bad pixels".format(cmd_args.exclude_regions_from_file)

# reset status from good -> bad in all specially earmarked regions
if cmd_args.include_regions_from_file:
    regions = pyregion.open(cmd_args.include_regions_from_file)
    mask = regions.get_mask(hdu=inhdu)
    badpix = badpix | mask
    if cmd_args.verbose:
        print "Regions from {} added to bad pixels".format(cmd_args.include_regions_from_file)

if cmd_args.verbose:
    nbad = badpix.sum()
    print "Number of bad pixels: {} ({:.3f}% of total)".format(nbad, float(100*nbad)/inhdu.data.size)


# Save the bad pixels candidates (some of these may be rejected later)
if cmd_args.debug or cmd_args.onlybadpix:
    save_aux_image(badpix.astype(int), "candidates")

if cmd_args.onlybadpix:
    sys.exit()                  # In this case, our work is done


## Now, label a list of contiguous regions of bad pixels (objects)
# Use 8 neighbours (including corners) instead of the default 4
labels, nlabels = ndi.label(badpix, structure=np.ones((3,3)))
if cmd_args.debug:
    save_aux_image(labels, "labels")
## Find slices corresponding to each labeled region
objects = ndi.find_objects(labels)
if cmd_args.verbose:
    print "Number of distinct bad pixel objects found:", len(objects)


## For each object, replace bad pixels with an average of nearby good ones
nskipped = 0
objecttable = list()
for i, (yslice, xslice) in enumerate(objects):
    label = i+1
    mx = xslice.stop - xslice.start
    my = yslice.stop - yslice.start
    # Object is deemed too large if both x and y sizes exceed dmax
    isTooBig = mx > cmd_args.dmax and my > cmd_args.dmax
    isReallyTooBig =(mx > cmd_args.ddmax and my > cmd_args.ddmax) or (mx*my > 2*cmd_args.ddmax**2)
    # Object is deemed compact if more than a certain fraction of the
    # pixels in the enclosing rect are bad.  For a circle, this would be
    # pi/4 = 0.785
    fillfactor = float(badpix[(yslice,xslice)].sum())/(mx*my)
    isCompact = fillfactor >= 0.8
    # Also require that box is not too elongated
    aspect = float(max(mx,my))/min(mx,my)
    isCompact = isCompact and aspect < 1.3
    isFilament = fillfactor <= 0.3 or max(mx,my) > 5*min(mx,my)

    # We may not have enough good pixels in the
    # enclosing rectangle, so expand it by one pixel in all 4
    # directions
    box = (
        slice(yslice.start-1, yslice.stop+1, yslice.step),
        slice(xslice.start-1, xslice.stop+1, xslice.step)
    )

    # Calculate the average value of "good" pixels in the expanded box
    replace_value = np.mean(inhdu.data[box][~badpix[box]])
    # Check if the bad pixels are darker on average than the surroundings
    contrast = np.mean(inhdu.data[box][badpix[box]]) / replace_value
    isShadow = contrast < 0.9

    isSkipped = (isShadow and (not cmd_args.allow_shadows)) or \
        (isTooBig and isCompact) \
        or isReallyTooBig \
        or (cmd_args.reject_filaments and isFilament)

    objecttable.append([label, mx, my, 
                        isTooBig, isCompact, isShadow, 
                        fillfactor, aspect, contrast, 
                        " *"[isSkipped]])

    if isSkipped:
        # Assume this object is a star or other "real" object
        # Undo all the bad pixels associated with this object
        badpix[box][labels[box] == label] = False
        # Then skip
        nskipped += 1
    else:
        # If this one is a keeper, then replace all the bad pixels in
        # box with the average non-bad value in same box
        inhdu.data[box] = np.where(badpix[box], replace_value, inhdu.data[box])


# Write out the final set of bad pixels
save_aux_image(badpix.astype(int), "badpix")

# Write out the table of object properties
with open(outprefix + "-objects.tab", "w") as f:
    f.write("\t".join(["# label", "mx", "my", 
                       "TooBig", "Compact", "Shadow", 
                       "fillfac", "aspect", "contrast", 
                       "skip?"]) + "\n")
    for row in objecttable:
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{}".format(*row) + "\n")

if cmd_args.verbose:
    print "Number of objects skipped: ", nskipped
    print "Replacement of bad pixels complete"

## Finally, write out the result to a new fits file
if cmd_args.multi_hdu:
    hdulist.writeto(outprefix + ".fits", clobber=True)
else:
    inhdu.writeto(outprefix + ".fits", clobber=True)
