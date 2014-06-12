cat $1 | grep Parameters --after 1 | grep '\[' | awk '{ print $2 }' | awk -F, '{ print $1 }' > Vmax

cat $1 | grep Parameters --after 1 | grep '\[' | awk '{ print $3 }' > Km

octave script.m > /dev/null
