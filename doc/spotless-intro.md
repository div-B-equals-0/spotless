<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#sec-1">1. Spotless</a>
<ul>
<li><a href="#sec-1-1">1.1. Compatibility and dependencies</a></li>
<li><a href="#sec-1-2">1.2. Installation</a></li>
<li><a href="#sec-1-3">1.3. Usage</a>
<ul>
<li><a href="#sec-1-3-1">1.3.1. File handling</a></li>
<li><a href="#sec-1-3-2">1.3.2. Specifying the pixel-rejection method to use</a></li>
<li><a href="#sec-1-3-3">1.3.3. Pruning the list of bad pixel objects</a></li>
<li><a href="#sec-1-3-4">1.3.4. Manually excluding and including regions</a></li>
<li><a href="#sec-1-3-5">1.3.5. Which images are saved</a></li>
<li><a href="#sec-1-3-6">1.3.6. Examples</a></li>
</ul>
</li>
<li><a href="#sec-1-4">1.4. Command line reference</a></li>
</ul>
</li>
</ul>
</div>
</div>

# Spotless<a id="sec-1" name="sec-1"></a>

A python program for removing cosmic rays, stitching artifacts, and (optionally) stars from an astronomical image or spectrum

The cosmic ray removal process has three stages:

1.  Pixel rejection and grouping contiguous bad pixel regions into "objects"
2.  Pruning the list of objects in order to minimise false positives
3.  "Zapping" all the pixels in each object in the final list by replacing bad pixels by the average of the "good" neighbors

## Compatibility and dependencies<a id="sec-1-1" name="sec-1-1"></a>

Spotless has been tested on Python 2.7, 3.3, and 3.4. 

It depends on the following third-party libraries:

-   [Numpy](http://www.numpy.org) and [Scipy](http://www.scipy.org/install.html)
    -   If you are using a full-featured scientific python distribution such as [Anaconda](http://continuum.io/downloads.html) then you already have these.
-   [Scikits Image](http://scikit-image.org)
    -   `conda install scikits-image`
-   Either [Astropy](http://www.astropy.org) (recommended) or PyFITS (no longer maintained)
    -   `conda install astropy`
-   [Pyregion](https://pypi.python.org/pypi/pyregion)
    -   `pip install pyregion`
    -   At time of writing (27 Nov 2014) this is not available in the Anaconda repository, but it can be installed from the Python Package Index using `pip` or `easy_install`

## Installation<a id="sec-1-2" name="sec-1-2"></a>

Spotless consists of a single executable python script `spotless`.   You can copy to somewhere on your `PATH` (or your working directory), or just leave it where it is, in which case you will have to specify the full path every time you run it.

## Usage<a id="sec-1-3" name="sec-1-3"></a>

Basic usage: 

    spotless INPUT_FILE

where `INPUT_FILE` is the name (without the `.fits` suffix) of the FITS file that you want to process.  An output file named `INPUT_FILE-cr.fits` will be produced.  

### File handling<a id="sec-1-3-1" name="sec-1-3-1"></a>

By default, spotless will process the HDU named `SCI` in the input
file if it is found, otherwise the first HDU will be used.  A
different HDU may be specified by using the `--hdu-name` option.

The extra string that is added to the name of the output file (default: `cr`) may be changed using the `--output-id` option.

### Specifying the pixel-rejection method to use<a id="sec-1-3-2" name="sec-1-3-2"></a>

Three different methods are available for choosing which pixels to reject.  In all cases, the result is a list of candidate "bad pixel objects", each of which is a contiguous region of pixels that are deemed to be bad.  This list is then pruned according to additional criteria discussed below to give a final list of bad pixel objects. 

1.  Simple threshold

    This method simply rejects all pixels with values greater than a certain value.  
    To use *only* this method, specify the options `--method thresh --threshold VALUE` where `VALUE` is the value of the threshold.   The option `--threshold VALUE` can also be used at the same time as another method, in which case the two methods are combined. 

2.  Edge detection followed by hole filling

    This is currently the default method so no option is required, but it can also be explicitly selected with `--method edge=`.  The option `--data-range MIN MAX` should be used to specify the approximate minimum and maximum values of the image brightness.  There is also an option `--use-log-scale`, which sometimes improves the results for images with a large dynamic range.  

3.  Segmentation via watershed transform of the elevation map

    This method does not seem to work very well, and so is not recommended.  If you want to try it anyway, the option `--method segment` should be specified. 

### Pruning the list of bad pixel objects<a id="sec-1-3-3" name="sec-1-3-3"></a>

The pixel-rejection methods all tend to produce a large number of false positives.  Many of these can be effectively pruned with a few simple criteria.  Each candidate object is checked to make sure it is sufficiently small, sufficiently circular, and not a dark region.  The details of these checks can be controlled with the options `--dmax`, `--ddmax`, `--reject-filaments`, and `--allow-shadows`.  See the Command line reference (See section 1.4) below for details. 

### Manually excluding and including regions<a id="sec-1-3-4" name="sec-1-3-4"></a>

The option `--exclude-regions-from-file FILENAME.reg` allows you to specify a DS9 region file for exclusion from the algorithm.  All regions in this file will be marked as "good".  That is, no cosmic ray detection will be performed inside the regions.  It is often necessary to use this option to protect known stars in your image (unless you want the stars to be zapped). 

Contrariwise, the option `--include-regions-from-file FILENAME.reg` allows you to specify known "bad" regions.  All regions in this file will be automatically added to the final list of bad pixel objects to be zapped.  

### Which images are saved<a id="sec-1-3-5" name="sec-1-3-5"></a>

By default, the `SCI` HDU in the output file has had all its final bad pixels zapped (but if the option `--only-bad-pix` is specified then it is just copied through unchanged from the input file).  A further HDU named `badpix` is  written to the output file, which is a logical mask of the pixels that were zapped. Any additional HDUs in the input file are copied through unchanged. 

If the option `--debug` is specified then several extra image HDUs are saved to the output file.  It the case of the `edge` method these are: `scaled`, `edges`, `candidates`, and `labels`. 

In addition, a file `INPUT_FILE-objects.tab` is written with a table that lists the label number (corresponding to the `labels` image), sizes, and pruning criteria for all the candidate bad pixel objects. 

### Examples<a id="sec-1-3-6" name="sec-1-3-6"></a>

    spotless  --data-range 0 10.0 --allow-shadows --output-id cr --verbose --debug F547M

## Command line reference<a id="sec-1-4" name="sec-1-4"></a>


    ./spotless --help

    usage: spotless [-h] [--hdu-name HDU_NAME] [--output-id OUTPUT_ID]
                    [--method {thresh,edge,segment}] [--onlybadpix]
                    [--threshold THRESHOLD] [--dmax DMAX] [--ddmax DDMAX]
                    [--data-range MIN MAX] [--use-log-scale]
                    [--segment-pars LO HI]
                    [--edge-pars SIGMA LOW_THRESHOLD HIGH_THRESHOLD]
                    [--thick-edges] [--reject-filaments] [--allow-shadows]
                    [--clip-negative]
                    [--exclude-regions-from-file EXCLUDE_REGIONS_FROM_FILE]
                    [--include-regions-from-file INCLUDE_REGIONS_FROM_FILE]
                    [--verbose] [--debug] [--multi-hdu]
                    fitsfile
    
    Remove cosmic rays and other bad pixels from an image
    
    positional arguments:
      fitsfile              Name of input image FITS file (sans extension)
    
    optional arguments:
      -h, --help            show this help message and exit
      --hdu-name HDU_NAME   Which HDU to use from the FITS file (default: SCI)
      --output-id OUTPUT_ID
                            Extra string to add to output filename to
                            differentiate from the input file (default: cr)
      --method {thresh,edge,segment}
                            Algorithm to use to find the bad pixels (default:
                            edge)
      --onlybadpix          Only calculate the bad pixel map - do not replace
                            pixels in the image (default: False)
      --threshold THRESHOLD
                            Assume any pixel above this level is bad (default:
                            None)
      --dmax DMAX           Maximum diameter of features to zap. Leave alone any
                            roughly circular objects that are larger than this.
                            (default: 5)
      --ddmax DDMAX         Absolute maximum diameter of features to zap. Leave
                            alone any objects that are larger than this, whatever
                            their shape may be. (default: 10)
      --data-range MIN MAX  Range for data scaling (default: None)
      --use-log-scale       Use logarithmic data scaling (default: False)
      --segment-pars LO HI  For 'segment' method only: thresholds of scaled data
                            to seed the good/bad regions (default: (0.0, 1.0))
      --edge-pars SIGMA LOW_THRESHOLD HIGH_THRESHOLD
                            For 'edge' method only: parameters for the Canny edge
                            detection algorithm. See: http://scikits-
                            image.org/docs/dev/auto_examples/plot_canny.html
                            (default: (1.0, 0.1, 0.2))
      --thick-edges         Make the edges be 3 pixels wide instead of the default
                            1 (default: False)
      --reject-filaments    Try to reject objects that look filamentary, since
                            they are probably not cosmic rays (default: False)
      --allow-shadows       Also remove objects that are darker than their
                            surroudings (not generally advised, especially if you
                            have dark globules in your image!) (default: False)
      --clip-negative       Also remove all negative pixels (default: False)
      --exclude-regions-from-file EXCLUDE_REGIONS_FROM_FILE
                            Read DS9 regions from a file, which are to be marked
                            as definite good pixels (default: None)
      --include-regions-from-file INCLUDE_REGIONS_FROM_FILE
                            Read DS9 regions from a file, which are to be marked
                            as definite bad pixels (default: None)
      --verbose, -v         Print informative progress messages (default: False)
      --debug, -d           Save auxiliary images of intermediate steps (default:
                            False)
      --multi-hdu, -m       Only provided for backward compatibility - this
                            behavior is now the default. Work in multi-HDU mode.
                            This assumes that the image is in the "SCI" HDU in the
                            input file (the argument --hdu-index is ignored). All
                            additional HDUs in the input file are copied through
                            to the output file. Only one output file is written,
                            all auxilliary arrays ("edges", "labels", "badpix",
                            etc) are written as additional HDUs in the same file.
                            (default: True)
