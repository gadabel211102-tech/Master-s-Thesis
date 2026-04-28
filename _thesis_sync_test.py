from zipfile import ZipFile, ZIP_DEFLATED
from xml.etree import ElementTree as ET
from pathlib import Path

src = Path('/home/gadeaalonsoj/tfm/thesis_sync_test.docx')
out = Path('/home/gadeaalonsoj/tfm/thesis_sync_test.updated.docx')
ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
ET.register_namespace('w', ns['w'])

replacements = {
"Annotated VCFs were then merged into a canonical project-level variant workbook. This unified analytical layer preserved sample-level molecular context while providing the common input used for descriptive review, common-SNP definition, mapping, and downstream statistical analyses. The refreshed BED-versus-detected comparison step inherited the same representative consequence, impact, and predictor fields so that design-matching reports and downstream interpretation tables remained synchronised.": "Annotated VCFs were then merged into a canonical project-level annotation workbook. This workbook preserves sample-linked molecular context for biological interpretation, but it should be understood as an annotation resource rather than as a complete sample-by-locus genotype matrix. It therefore provides the common input used for descriptive review, common-SNP definition, mapping labels, and downstream interpretation of gene, consequence, impact, and predictor fields. The refreshed BED-versus-detected comparison step inherits the same representative consequence, impact, and predictor fields so that design-matching reports and interpretation tables remain synchronised, while explicit genotype-state logic is sourced separately from the forced-genotype matrix.",
"Descriptive analyses summarised the variant landscape across the targeted GSDMB locus and the surrounding chromosome 17 panel, including consequence and impact distributions, gene-level mapping summaries, and GSDMB-focused landscape plots. For the main inferential SNP branch, common variants were defined using a combined non-Finnish European gnomAD allele-frequency field, in which gnomADe_NFE_AF was used preferentially and gnomADg_NFE_AF was used as fallback when necessary. Variants were retained for the common-SNP backbone when the NFE allele frequency exceeded 0.01.": "Descriptive analyses summarised the variant landscape across the targeted GSDMB locus and the surrounding chromosome 17 panel, including consequence and impact distributions, gene-level mapping summaries, and GSDMB-focused landscape plots. Annotation labels for these summaries were taken from the annotation workbook, whereas genotype-state plots used the explicit forced-genotype matrix so that wild-type, heterozygous ALT, homozygous ALT, and uncallable states were kept analytically distinct. For the main inferential SNP branch, common variants were defined using a combined non-Finnish European gnomAD allele-frequency field, in which gnomADe_NFE_AF was used preferentially and gnomADg_NFE_AF was used as fallback when necessary. Variants were retained for the common-SNP backbone when the NFE allele frequency exceeded 0.01.",
"SNP-clinical association analyses were performed in tumour samples only, with breast and endometrial cohorts analysed separately because the available clinical variables differed between the two datasets. Continuous outcomes were assessed with Mann-Whitney U tests, categorical outcomes with Fisher's exact or chi-square tests as appropriate, and genotype-aware summaries further examined wild-type, heterozygous, and homozygous ALT states where sample support was sufficient. The refreshed genotype logic used the same ALT-state coding across descriptive summaries, focused genotype matrices, and inferential analyses so that follow-up interpretations remained directly comparable across workflow stages.": "SNP-clinical association analyses were performed in tumour samples only, with breast and endometrial cohorts analysed separately because the available clinical variables differed between the two datasets. Continuous outcomes were assessed with Mann-Whitney U tests, categorical outcomes with Fisher's exact or chi-square tests as appropriate, and genotype-aware summaries further examined wild-type, heterozygous, and homozygous ALT states where sample support was sufficient. In the refreshed workflow, these genotype states were taken from the explicit forced-genotype matrix rather than reconstructed from annotation-derived carrier counts. Callable percentages therefore used locus-specific callable denominators, while low-depth and no-call samples were kept separate instead of being absorbed into the wild-type group. The same ALT-state coding was then carried across descriptive summaries, focused genotype matrices, and inferential analyses so that follow-up interpretations remained directly comparable across workflow stages.",
}

with ZipFile(src) as zin:
    root = ET.fromstring(zin.read('word/document.xml'))
    changed = 0
    for p in root.findall('.//w:p', ns):
        ts = p.findall('.//w:t', ns)
        full = ''.join([t.text for t in ts if t.text])
        if full in replacements:
            ts[0].text = replacements[full]
            for t in ts[1:]:
                t.text = ''
            changed += 1
    print('changed', changed)
    new_xml = ET.tostring(root, encoding='utf-8', xml_declaration=True)
    with ZipFile(out, 'w', ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = new_xml if item.filename == 'word/document.xml' else zin.read(item.filename)
            zout.writestr(item, data)

with ZipFile(out) as z:
    root = ET.fromstring(z.read('word/document.xml'))
paras = []
for p in root.findall('.//w:p', ns):
    texts = [t.text for t in p.findall('.//w:t', ns) if t.text]
    if texts:
        para = ''.join(texts).strip()
        if para:
            paras.append(para)
for i in [251, 256, 261]:
    print(f'[{i}] {paras[i]}')
