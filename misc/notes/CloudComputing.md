# Cloud computing
Here I provide a guide for users that want to use the pipeline to analyse public SRA/GEO data or their own fastq files. This is under construction and will undergo major changes.

## Google cloud
Open a paid Google Cloud account - as a news users receive a free credit of $300. At the dashboard, click the menu button on the top left of screen and click "Compute Engine". You might have to wait a few minutes. Click "Create an instance template" and follow the process to specify the the image characteristics suitable for mamalian sized genomes:
* 16+ vCPUs
* 60 GB memory
* OS:Default Debian is OK or Ubuntu 16.04
* 60 to 100 GB disk
* Select preemtibility - this will save you ~70% of the compute cost but may be intermittent service
* Once the vm is launched, on the dahboard you can click "SSH". From here you can follow the "quick start guide" I posted on the main page README.md
### Optional steps for repeat/power users
* Specify startup command. Here is an example. This may require a bit of testing to get working.
```while true ; do docker run mziemann/tallyup mmusculus >dee2.log 2>&1``` 
* Save the image template to enable future launching. Save the command line (very bottom of page) to enable CLI interface using the "gcloud" command.

## Amazon AWS
Launch an EC2 instance. Give it a name. Select Ubuntu or Amazon Linux. Choose an instance type based on the memory requirements, which are linked to the reference genome size. Initially, try m4.xlarge for testing. Next click "configure instance details". Click spot price and set maximum price. Continue to disk, attach 60 GB of storage. Security group? Launch instance. Proceed without a key pair, you need to remember the account password entered previously. note the address to enable ssh login the 
