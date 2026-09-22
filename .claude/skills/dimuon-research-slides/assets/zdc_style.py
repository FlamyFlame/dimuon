"""
zdc_style.py — shared slide style kit for Yuhan's interview deck.

Import this in any deck-building script so every section (ZDC calibration,
ZDC online monitoring, dimuon analysis) shares ONE look: same palette,
fonts, header treatment, cards, stat callouts, and figure handling.

Usage:
    from zdc_style import *
    prs = new_deck()
    s = add_slide(prs); bg(s, WHITE)
    header(s, "The task", "My slide title")
    bullets(s, 0.7, 1.95, 12.4, 2.0, [
        [["Plain text, ", DKGRAY], ["highlighted", BLUE, True], [".", DKGRAY]],
        "A simple one-color bullet.",
    ], size=17)
    pic_fit(s, "figs/plot.png", 0, 3.9, 12.0, 3.3, center_x=SW/2)
    finish(prs, "MyDeck.pptx")          # saves + rezips

Conventions that keep sections consistent (follow these):
  * Layout = 13.333 x 7.5 in (16:9 wide). SW, SH below.
  * Title slide on NAVY; content slides on WHITE; optional section divider on NAVY.
  * Slide titles are ONE line (34pt). If a title wraps, shorten it.
  * Fill the page: prefer full-width figures along the bottom over a small
    figure floating in a half-empty column.
  * Body 16-18pt, captions 13pt italic muted, stat numbers 44pt.
  * Motif: three accent dots (steel/red/green) on the title; rounded cards
    with a light tint + soft shadow. NEVER use accent bars/stripes/underlines.
  * Always add speaker notes (s.notes_slide.notes_text_frame.text = "...").
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

# ------------------------------------------------------------------ palette
NAVY   = RGBColor(0x1E, 0x2A, 0x44)   # title/divider background, slide titles
BLUE   = RGBColor(0x2E, 0x5A, 0x88)   # primary highlight
STEEL  = RGBColor(0x4C, 0x78, 0xA8)   # kicker / bullet dots / accent
RED    = RGBColor(0xC0, 0x39, 0x2B)   # problem / "the BRAN" / negative
GREEN  = RGBColor(0x2E, 0x7D, 0x32)   # fixed / positive / confirmation
GRAY   = RGBColor(0x5B, 0x66, 0x75)   # muted captions/labels
DKGRAY = RGBColor(0x37, 0x41, 0x4F)   # body text
LIGHT  = RGBColor(0xF2, 0xF5, 0xF9)   # card tint
WHITE  = RGBColor(0xFF, 0xFF, 0xFF)
ICE    = RGBColor(0xCA, 0xDC, 0xFC)   # subtitle on navy

HEAD_FONT = "Cambria"     # titles / formulas / stat numbers (safe-list serif)
BODY_FONT = "Calibri"     # body / bullets / captions (safe-list sans)

SW, SH = 13.333, 7.5

# ------------------------------------------------------------------ core
def new_deck(author="Yuhan Guo", title="ZDC / dimuon"):
    prs = Presentation()
    prs.slide_width  = Inches(SW)
    prs.slide_height = Inches(SH)
    prs.author = author
    prs.title = title
    return prs

def add_slide(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])  # blank

def bg(s, color):
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = color

def box(s, x, y, w, h):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    return tb, tf

def setrun(p, text, size, color, bold=False, italic=False, font=BODY_FONT):
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = bold; r.font.italic = italic
    r.font.name = font; r.font.color.rgb = color
    return r

def notes(s, text):
    s.notes_slide.notes_text_frame.text = text

# ------------------------------------------------------------------ blocks
def header(s, kicker, title, kcolor=STEEL, tcolor=NAVY):
    """Kicker (uppercase, spaced) + one-line title. No underline rule."""
    _, tf = box(s, 0.7, 0.42, 12.2, 0.45)
    p = tf.paragraphs[0]
    setrun(p, kicker.upper(), 15, kcolor, bold=True, font=BODY_FONT)
    p.runs[0].font._rPr.set('spc', '200')
    _, tf2 = box(s, 0.7, 0.82, 12.4, 1.0)
    setrun(tf2.paragraphs[0], title, 34, tcolor, bold=True, font=HEAD_FONT)

def bullets(s, x, y, w, h, items, size=17, color=DKGRAY, gap=9):
    """items: each item is either a string, or a list of [text, color, bold?] segments."""
    _, tf = box(s, x, y, w, h)
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(gap); p.space_before = Pt(0)
        setrun(p, "•  ", size, STEEL, bold=True)
        if isinstance(it, list):
            for seg in it:
                setrun(p, seg[0], size, seg[1] if len(seg) > 1 else color,
                       bold=(len(seg) > 2 and seg[2]))
        else:
            setrun(p, it, size, color)

def caption(s, x, y, w, text, size=13):
    _, tf = box(s, x, y, w, 0.6)
    setrun(tf.paragraphs[0], text, size, GRAY, italic=True)

def card(s, x, y, w, h, fill=LIGHT, shadow=True):
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                            Inches(x), Inches(y), Inches(w), Inches(h))
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    sh.line.fill.background()
    sh.shadow.inherit = False
    if shadow:
        el = sh.shadow._element
        # subtle outer shadow
        from pptx.oxml.ns import qn
        sp = sh._element.spPr
        efflst = sp.makeelement(qn('a:effectLst'), {})
        sh_el = sp.makeelement(qn('a:outerShdw'),
                               {'blurRad': '60000', 'dist': '20000', 'dir': '5400000'})
        clr = sp.makeelement(qn('a:srgbClr'), {'val': '000000'})
        alpha = sp.makeelement(qn('a:alpha'), {'val': '12000'})
        clr.append(alpha); sh_el.append(clr); efflst.append(sh_el); sp.append(efflst)
    return sh

def stat(s, x, y, w, big, label_lines, bigcolor, h=1.7):
    """Big number callout with up to 2 centered label lines."""
    card(s, x, y, w, h)
    _, tf = box(s, x+0.1, y+0.16, w-0.2, 0.85); tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    setrun(p, big, 44, bigcolor, bold=True, font=HEAD_FONT)
    _, tf2 = box(s, x+0.08, y+1.04, w-0.16, 0.6)
    for i, ln in enumerate(label_lines):
        p2 = tf2.paragraphs[0] if i == 0 else tf2.add_paragraph()
        p2.alignment = PP_ALIGN.CENTER; p2.space_after = Pt(0); p2.space_before = Pt(0)
        setrun(p2, ln, 13, GRAY)

def pic_fit(s, path, x, y, maxw, maxh, center_x=None):
    """Place an image preserving aspect ratio, capped to maxw x maxh.
    Returns (x, y, w, h). Use center_x=SW/2 to center horizontally."""
    iw, ih = Image.open(path).size
    ar = iw / ih
    w = maxw; h = w / ar
    if h > maxh:
        h = maxh; w = h * ar
    if center_x is not None:
        x = center_x - w / 2
    s.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))
    return x, y, w, h

# ------------------------------------------------------------------ slide presets
def title_slide(prs, title, subtitle, footer):
    s = add_slide(prs); bg(s, NAVY)
    _, tf = box(s, 0.9, 1.95, 11.7, 2.6)
    setrun(tf.paragraphs[0], title, 50, WHITE, bold=True, font=HEAD_FONT)
    p2 = tf.add_paragraph(); p2.space_before = Pt(16)
    setrun(p2, subtitle, 23, ICE)
    for i, c in enumerate([STEEL, RED, GREEN]):
        o = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.95+i*0.36), Inches(4.75),
                               Inches(0.22), Inches(0.22))
        o.fill.solid(); o.fill.fore_color.rgb = c
        o.line.fill.background(); o.shadow.inherit = False
    _, tf3 = box(s, 0.9, 6.35, 11.7, 0.6)
    setrun(tf3.paragraphs[0], footer, 16, RGBColor(0x9F, 0xB0, 0xC8))
    return s

def section_divider(prs, kicker, title):
    s = add_slide(prs); bg(s, NAVY)
    _, tf = box(s, 0.9, 2.7, 11.5, 0.5)
    setrun(tf.paragraphs[0], kicker.upper(), 15, STEEL, bold=True)
    tf.paragraphs[0].runs[0].font._rPr.set('spc', '200')
    _, tf2 = box(s, 0.9, 3.15, 11.5, 1.4)
    setrun(tf2.paragraphs[0], title, 40, WHITE, bold=True, font=HEAD_FONT)
    return s

def finish(prs, out_path):
    prs.save(out_path)
    print("saved", out_path)
