
FIGS="standard-model stock-spec-catch-model"

for i in $FIGS; do
  pdflatex $i.tex
  pdfcrop $i.pdf
  magick -density 600 $i-crop.pdf -background white -alpha remove -alpha off -quality 100 $i-crop.png
  cp $i-crop.pdf $i-crop.png ../images/
done
