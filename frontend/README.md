# Frontend scripts
This folder describes scripts and files that are required for the operation of the public facing web server. 

## cgibin
This folder of scripts is to be placed in the /usr/lib/cgi-bin directory.
The <b>acc.sh</b> script manages the queue of accessions that are requested by the worker nodes. 
The <b>search.sh</b> script as the name suggests, searches metadata for matching DEE2 datasets. 
The <b>request.sh</b> script packages the datasets for the end user and bundles it up in a zip archive.
The <b>request2.sh</b> script sends the STAR gene count matrix to [Degust](http://degust.erc.monash.edu) to enable differential analysis in the browser.
Similarly <b>request3.sh</b> sends the Kallisto transcript level counts to Degust.
Lastly <b>request4.sh</b> sends the Kallisto counts sumarised to gene level to Degust.

# html
This folder contains html files that users interact with through the browser.

# clean.sh
This is the script that performs initial checking of incoming datasets and relays them to the data repo.
