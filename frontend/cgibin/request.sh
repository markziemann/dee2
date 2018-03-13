#!/bin/bash
QUERY_STRING2=$(echo $QUERY_STRING | sed 's/\&/\\&/g')
#echo $QUERY_STRING2
ssh -i ~/.ssh/monash/id_rsa mziemann@118.138.246.227 ./extract.sh $QUERY_STRING2

