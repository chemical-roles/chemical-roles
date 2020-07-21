![Chemical Roles Logotype](CRoG-logotype-1024.png)

# The Chemical Roles Logo

## Design

The Chemical Roles logo consists of a stylized frog's foot overlaid with a
representation of a graph and the text "CRoG," which stands for "Chemical ROles
Graph." The use of the frog's foot takes inspiration from the pronunciation of
"CRoG" sounding like something a frog might say.

The text is set in Trattatello Regular with manual kerning and adjustments.

The alternative "CRoG King" logo includes the same elements as the main logo
superimposed on the icon of a frog wearing a crown.

## Files

The CRoG logos are provided in several formats and sizes:

- `CRoG-logotype.svg`, an SVG of the logo with the text stored as a `<text>`
  element
- `CRoG-logotype-paths.svg`, an SVG of the logo with the text converted to
  paths, allowing the logo to be displayed on systems wich don't have the
  necessary font
- `CRoG-square.svg`, an SVG of the logo without the text
- `CRoG-{logotype, square}-{100, 1024}.png`, rasterized versions of the logotype
  and logo at two sizes

There are similar versions for the CRoG king logo, in both a vertical and
horizontal format.

## Creating Graphics

The CRoG logos were created as vector graphics in Inkscape and rasterized to PNG
with Inkscape as well. The rasterized `*-square` versions were made actually
square by adding transparent padding with this ImageMagick command
([reference](https://imagemagick.org/Usage/thumbnails/#pad)):

```sh
for f in *square*.png
do
    magick "$f" -virtual-pixel none \
        -set option:distort:viewport "%[fx:max(w, h)]x%[fx:max(w, h)]-%[fx:max((h-w)/2, 0)]-%[fx:max((w-h)/2, 0)]" \
        -filter point -distort ScaleRotateTranslate 0 +repage "$f"
done
```

The PNG versions were then [optimized](https://blog.codinghorror.com/zopfli-optimization-literally-free-bandwidth/)
using `optipng` and `advdef`:

```sh
optipng -o2 -nb *.png && advdef -z -4 *.png
```

## License

The CRoG logos are derivatives of
"[Frog foot icon](https://game-icons.net/1x1/delapouite/frog-foot.html)"
by [Delapouite](http://delapouite.com/), used under
[CC BY 3.0](http://creativecommons.org/licenses/by/3.0/).

The CRoG King logos are derivatives of the "Frog foot icon" by Delapouite as
above, as well as
"[Frog prince icon](https://game-icons.net/1x1/delapouite/frog-prince.html)"
also by Delapouite, used under
[CC BY 3.0](http://creativecommons.org/licenses/by/3.0/).

The logos are licensed under
[CC BY 3.0](http://creativecommons.org/licenses/by/3.0/)
by [Scott Colby](https://github.com/scolby33).

<p align="center">
    <img src="CRoG-king-100.png" alt="The CRoG King">
</p>
