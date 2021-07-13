#!/bin/bash

ls *aa.tab > tmp.tmp
echo "#/!bin/bash">tmp2.sh

while read -r ONE; do

    POS=$(echo "$ONE" | cut -d_ -f2)
    #echo $POS
    OFFSET=241
    OPOS=$((POS - OFFSET))
    echo "sed -n \"$OPOS p\" $ONE" >> tmp2.sh

done < tmp.tmp
chmod 755 tmp2.sh
tmp2.sh
rm tmp.tmp
#rm tmp2.sh