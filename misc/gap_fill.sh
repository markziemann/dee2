#!/bin/bash

for ORG in athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae ; do
  MDATA=../sradb/${ORG}_metadata.tsv.cut
  DATA=../data/${ORG}
  for SRP in $(sed 1d $MDATA | cut -f5 | sort -u  ) ; do
    TOT=0 ; DONE=0 ; NOT_DONE=""
    for SRR in $(grep -w $SRP $MDATA | cut -f1 ) ;  do
      TOT=$((TOT+1))
      if [ -d $DATA/$SRR ] ; then
        DONE=$((DONE+1))
      else
        NOT_DONE=${NOT_DONE},$SRR
      fi
    done
    echo $SRP $TOT $DONE $NOT_DONE
  done | awk '$3/$2<1 && $3/$2>.5 {print $1,$2,$3,$3/$2,$4}' | sort -k5g \
  |cut -d ' ' -f5 | tr ',' '\n' | sed '/^$/d' > ${ORG}.queue.txt

  LENGTH=$(wc -l <  ${ORG}.queue.txt )
  if [ $LENGTH -gt 4000 ] ;  then
    awk '{print $0, $0}' ${ORG}.queue.txt | cut -c4- | sort -k1gr | awk '{print $2}' > tmp
    mv tmp ${ORG}.queue.txt
    split -n4 ${ORG}.queue.txt ${ORG}.queue.txt.
    scp ${ORG}.queue.txt.aa mdz@10.115.170.32:/home/mdz/bfx/dee2/gap_fill
    scp -i ~/.ssh/monash/id_rsa ${ORG}.queue.txt.ab mziemann@118.138.234.73:/home/mziemann/dee2/gap_fill/
    scp -i ~/.ssh/monash/id_rsa ${ORG}.queue.txt.ac mziemann@118.138.234.56:/home/mziemann/dee2/gap_fill/
    scp -i ~/.ssh/monash/id_rsa ${ORG}.queue.txt.ad mziemann@118.138.246.227:/scratch/mziemann/dee2/gap_fill/
  fi
done

exit
###
# Run this on the other machine
#!/bin/bash
ORG=$(echo $0 | sed 's@./gap_fill_@@' | sed 's@.sh$@@')
ABBREV=$(echo $ORG | cut -c-3)
#check path to queue
for FILE in /scratch/mziemann/dee2/gap_fill/${ORG}.queue.txt.* ; do
  for SRR in $(shuf $FILE ) ; do
    echo $SRR
    docker run mziemann/tallyup_${ABBREV} $ORG $SRR && sed -i "/^${SRR}/d" $FILE
  done
done

