# DEE2
The aim of DEE2 is to make all RNA-seq data freely available to everyone. DEE2 consists of three parts:
* Webserver where end-users can search for and obtain data-sets of interest 
* Pipeline that can download and process SRA data as well as users' own fastq files.
* Back-end that collects, filters and organises data provided by contributing worker nodes.

DEE2 currently supports analysis of several major species including A. thaliana , C. elegans, D. melanogaster, D. rerio, E. coli, H. sapiens, M. musculus, R. norvegicus and S. cerevisiae. The DEE2 pipeline downloads data from SRA and processes it, providing tabulated data that can be used in downstream statistical analysis.

## How can I access the processed data?
The processed data is available at http://dee2.io and can be also accessed using our specially developed [R interface](../master/AccessDEEfromR.md)

***Due to planned maintenance of the Nectar research cloud, there will be an 8 hour outage some time
between 10am 19th June and 6pm 4th June (AEST).***

## Want to learn more?
For information on different parts of the app, see the specific documentation:
* [How to use the pipeline to process raw data](../master/pipeline/README.md)
* [How we validated the accuracy of the pipeline](../master/validation/README.md)
* [How the backend scripts are organised](../master/backend/README.md)
* [How the frontend scripts are organised](../master/frontend/README.md)

Feedback, bug reports, and contributions to code development are very welcome. 
