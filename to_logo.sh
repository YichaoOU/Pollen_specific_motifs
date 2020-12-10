module load meme

module load texlive/20190410

meme2images -eps all_15_motifs.pwm Motif_logo

cd Motif_logo

for i in *.eps; do sed -i 's/showFineprint true/showFineprint false/' $i;done

for i in *.eps;do epstopdf $i -o=$i.pdf;done


