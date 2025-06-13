# DEE2
The aim of DEE2 is to make all RNA-seq data freely available to everyone. DEE2 consists of three parts:
* Webserver where end-users can search for and obtain data-sets of interest 
* Pipeline that can download and process SRA data as well as users' own fastq files.
* Back-end that collects, filters and organises data provided by contributing worker nodes.

DEE2 currently supports analysis of several major species including A. thaliana , C. elegans, D. melanogaster, D. rerio, E. coli, H. sapiens, M. musculus, R. norvegicus and S. cerevisiae. The DEE2 pipeline downloads data from SRA and processes it, providing tabulated data that can be used in downstream statistical analysis.

We have recently added support for B. distachyon, G. max, H. vulgare, O. sativa, P. trichocarpa, S. bicolor, S. lycopersicum, S. tuberosum, T. aestivum, V. vinifera and Z. mays.

## How can I access the processed data?
The processed data is available at http://dee2.io and can be also accessed using our specially developed [R interface](https://bioconductor.org/packages/getDEE2/).
Project bundles are available [here](https://dee2.io/huge/), and bulk data dumps are available [here](https://dee2.io/mx/).

<s>If there is a particular dataset of interest missing from DEE2, you can use the [request webform](http://dee2.io/request.html) to have it expedited.</s>

## Want to learn more?
For information on different parts of the app, see the specific documentation:
* [How to use the pipeline to process raw data](../master/pipeline/README.md)
* [How we validated the accuracy of the pipeline](../master/validation/README.md)
* [How the backend scripts are organised](../master/backend/README.md)
* [How the frontend scripts are organised](../master/frontend/README.md)

Feedback, bug reports, and contributions to code development are very welcome. 
