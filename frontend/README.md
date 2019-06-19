# Frontend scripts
This folder describes scripts and files that are required for the operation of the public facing web server. 

## cgibin
This folder of scripts is to be placed in the /usr/lib/cgi-bin directory. The acc.sh script manages the queue of accessions that are requested by the worker nodes. The search.sh script as the name suggests, searches metadata for matching DEE2 datasets. The request.sh script packages the datasets for the end user and bundles it up in a zip archive. The request2.sh script sends the STAR gene count matrix to [Degust](http://degust.erc.monash.edu) to enable differential analysis in the browser.

# html
This folder contains html files that users interact with through the browser.

# clean.sh
This is the script that performs initial checking of incoming datasets and relays them to the data repo.
