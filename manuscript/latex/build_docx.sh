#!/usr/bin/env bash
# Build main.docx and supplementary.docx from the LaTeX sources with pandoc, and copy them
# into docs/ beside the PDFs. CI builds only the PDFs; run this by hand after latexmk.
#
# Word cannot display PDF figures, so the PNG copies from results/figures/ are swapped in.
# Citations are rendered by pandoc's citeproc (author-date), not by sn-nature.bst.
set -euo pipefail
cd "$(dirname "$0")"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures"
cp ../../results/figures/*.png ../../results/figures/supplementary/*.png "$tmp/figures/"
cp references.bib "$tmp/"

for f in main supplementary; do
    sed -E 's/(\\includegraphics(\[[^]]*\])?\{[^}]*)\.pdf\}/\1.png}/' "$f.tex" > "$tmp/$f.tex"
    if [ "$f" = supplementary ]; then
        # The SI points into the manuscript through xr, which pandoc cannot follow, so those
        # \ref{}s are resolved from main.aux here. latexmk must have built main.pdf first.
        perl -e 'open A, "main.aux" or die "main.aux missing: run latexmk on main.tex first\n";
                 while (<A>) { $n{$1} = $2 if /\\newlabel\{([^}]+)\}\{\{([^}]*)\}/ }
                 local $/; $_ = <STDIN>;
                 s/\\ref\{([^}]+)\}/exists $n{$1} ? $n{$1} : "\\ref{$1}"/ge; print' \
            < "$tmp/$f.tex" > "$tmp/$f.resolved" && mv "$tmp/$f.resolved" "$tmp/$f.tex"
    fi
    (cd "$tmp" && pandoc "$f.tex" --citeproc --bibliography=references.bib \
        --resource-path=.:figures -o "$f.docx")
    cp "$tmp/$f.docx" .
done

cp main.docx ../../docs/Photorespiration_multiomics_microgravity.docx
cp supplementary.docx ../../docs/Photorespiration_multiomics_microgravity_SI.docx
echo "wrote main.docx, supplementary.docx and their docs/ copies"
