import sys, os, math
import gelscore, translate
import numpy as np
import svgwrite, collections, colorsys
from svgwrite import cm, mm, rgb, deg

'''
Start a new drawing
Given a sequence and domain boundaries
	- Draw line for sequence
	- Draw box for domains
'''

def proteinName(x):
	return x[0].upper() + x[1:].lower()

def hexToRgb(x):
	return tuple(map(ord, x.decode('hex')))

def rgbToHex(rgbtuple):
	return "".join(map(chr, rgbtuple)).encode('hex')

def ycbcrToRgb(y, cb, cr):
    r = int(y + 1.402 * (cr-128))
    g = int(y - .34414 * (cb-128) -  .71414 * (cr-128))
    b = int(y + 1.772 * (cb-128))
    return (r, g, b)

def rgbToYcbcr(r, g, b):
    y = int(.299*r + .587*g + .114*b)
    cb = int(128 -.168736*r -.331364*g + .5*b)
    cr = int(128 +.5*r - .418688*g - .081312*b)
    return (y, cb, cr)

def princompScore(x, mat):
	return x.dot(mat)

def normVector(x):
	return x/np.sqrt((x*x).sum())

def v01To8bit(x):
	"""Convert floating-point vector on (0,1) to integer vector on (0,255)"""
	return tuple(np.int_(np.floor(255*np.r_[x])))

def vectorToColor(x, comps, colormap='rgb', whichcomps=slice(0,3)):
	"""Projects normalized vector x onto principal components comps"""
	pc = princompScore(x, comps)
	npc = normVector(pc[whichcomps])
	npc = (npc+1)/2.0 # Map to 0,1
	#print npc
	if colormap=='rgb':
		rgb = npc
	elif colormap=='hsv':
		rgb = colorsys.hsv_to_rgb(*npc)
		#print rgb
	elif colormap=='hls':
		rgb = colorsys.hls_to_rgb(*npc)
	elif colormap=='ycbcr':
		rgb = ycbcrToRgb(*v01To8bit(npc))
	rgb8bit = v01To8bit(rgb)
	#print rgb8bit
	res = rgbToHex(rgb8bit)
	return res

def sequenceCompositionNorm(seq):
	c = collections.Counter(seq)
	v = [0]*20
	for (i,aa) in enumerate(translate.AAs()):
		v[i] = c[aa]
	vn = normVector(np.array(v))
	return vn

def sequenceToColor(seq, comps, colormap='rgb', whichcomps=slice(0,3)):
	return vectorToColor(sequenceCompositionNorm(seq), comps, colormap, whichcomps)
	

class SequenceDrawing(object):
	cm_per_point = 0.035278 # cm per PostScript/font point
	
	def __init__(self, width, height, cm_per_aa=0.25, startx=2, starty=2, label_startx=1.0, domain_height=0.25, row_height=2.0):
		self._drawing = svgwrite.Drawing(size=(width*cm, height*cm))
		# Initialize scripts
		scripts = self._drawing.script(href="../src/domain-display.js")
		self._drawing.add(scripts)
		self._cm_per_aa = cm_per_aa
		self._startx = startx
		self._starty = starty
		self._label_startx = label_startx
		self._domain_height = domain_height
		self._row_height = row_height
		self._row = 0
	
	def write(self, fname):
		self._drawing.saveas(fname)
	
	def nextRow(self):
		self._row += 1
	
	@property
	def row(self):
		return self._row

	@row.setter
	def row(self, therow):
		self._row = therow
	
	@property
	def rowpos(self):
		return self._starty + self._row*self._row_height
	
	def cartoon(self, seq, gene_name, domain_boundaries=[], domain_colors=None, text_color='black'):
		dg = self._drawing
		self._seq = seq
		self._name = gene_name
		if domain_colors is None:
			domain_colors = ['gray']*len(domain_boundaries)
		domain_boundaries.sort()
		gp = dg.add(dg.g()) #stroke='black', stroke_linecap='round', stroke_width=1.0))
		self.lines(domain_boundaries, gp)
		self.boxes(domain_boundaries, domain_colors, gp)
		self.coordinates(domain_boundaries, gp)
		self.subsequences(domain_boundaries, gp)
	
	def label(self, labeltext, font='Arial', size=12, color='black', weight='normal'):
		self.text(labeltext, x=self._label_startx, y=self.rowpos + SequenceDrawing.cm_per_point*size/4.0, 
			size=size, font='Arial', color=color, weight=weight, align='right')
		
	def domainBounds(self, aastart, aaend):
		"""Return x1,y1, width, height"""
		cmperaa = self._cm_per_aa
		x = self._startx + cmperaa*(aastart-1)
		y = self.rowpos - self._domain_height/2.0
		width = cmperaa*(aaend - aastart)
		height = self._domain_height
		return (x, y, width, height)
	
	def lineBounds(self, aastart, aaend):
		"""Return x1, x2"""
		cmperaa = self._cm_per_aa
		x1 = self._startx + cmperaa*(aastart-1)
		x2 = x1 + cmperaa*(aaend-aastart)
		y1 = y2 = self.rowpos
		res = (x1, y1, x2, y2)
		#print aastart, aaend, res
		return res
	
	def lines(self, domain_boundaries, group):
		dg = self._drawing
		linegroup = dg.g(stroke_width=1.0, stroke='black', stroke_linecap='round')
		group.add(linegroup)
		if len(domain_boundaries) > 0:
			(x1s, y1s, x2s, y2s) = self.lineBounds(1, domain_boundaries[0][0])
			#print "\t",x1s, x2s
			for (st, en) in domain_boundaries:
				#print "\t",st, en
				(x1, y1, x2, y2) = self.lineBounds(st, en)
				#print "\t",x1, x2
				#print "\t", x1s, x1
				ln = dg.line(start=(x1s*cm, y1s*cm), end=(x1*cm, y1*cm))
				linegroup.add(ln)
				(x1s, y1s) = (x2, y2)
			# Finish up
			if en < len(self._seq):
				(x1, y1, x2, y2) = self.lineBounds(en, len(self._seq))
				ln = dg.line(start=(x1*cm, y1*cm), end=(x2*cm, y2*cm))
				linegroup.add(ln)
		else:
			(x1, y1, x2, y2) = self.lineBounds(1, len(self._seq))
			ln = dg.line(start=(x1*cm, y1*cm), end=(x2*cm, y2*cm))
			linegroup.add(ln)
			

	
	def boxes(self, domain_boundaries, domain_colors, group):
		dg = self._drawing
		# Make boxes
		boxgroup = dg.g(
			onmouseover="borderHighlight(evt); show(evt, evt.target.id, ['coords','seq'])", 
			onmousedown="toggleVisibility(evt);", 
			onmouseout ="borderRestore(evt); hide(evt, evt.target.id, ['coords','seq'])")
		group.add(boxgroup)
		for (i, (st, en)) in enumerate(domain_boundaries):
			col = domain_colors[i]
			self.box(st, en, col, boxgroup)
	
	def box(self, aastart, aaend, color, group=None):
		dg = self._drawing
		(x, y, width, height) = self.domainBounds(aastart, aaend)
		roundradius = self._domain_height
		subseq = "{:s}-{:d}-{:d}".format(self._name, aastart, aaend)
		domain_rect = dg.rect(
			id=subseq,
			insert=(x*cm, y*cm), size=(width*cm, height*cm), 
			rx = roundradius*mm, ry=roundradius*mm,	fill=color)
		if not group is None:
			group.add(domain_rect)
		else:
			dg.add(domain_rect)
	
	def text(self, text, x, y, group=None, font='Arial', size=12, color='black', weight='normal', align=None, display='inline'):
		dg = self._drawing
		align_dict = {None:None, 'center':'middle', 'left':'start', 'right':'end'}
		# DAD: don't write out styles if you don't need them
		#style_tag = "font-size:{size}px; font-family:{font}; font-weight:{weight}; \
		#	fill:{color}; text-anchor:{align}; alignment-baseline:{valign};".format(
		#	size=size, font=font, weight=weight, color=color, align=align_dict[align], valign=vertical_align_dict[valign])
		textgroup = dg.g(display=display)		
		textgroup.add(dg.text(text, (x*cm,y*cm), font_family=font, font_size=size, font_weight=weight, text_anchor=align_dict[align]))
		if group is None:
			dg.add(textgroup)
		else:
			group.add(textgroup)
	
	def coordinates(self, domain_boundaries, group=None, size=9):
		dg = self._drawing
		if group is None:
			group = dg		
		if len(domain_boundaries) > 0:
			for (st, en) in domain_boundaries:
				(x1, y1, x2, y2) = self.lineBounds(st, en)
				subseq = "{:s}-{:d}-{:d}-coords".format(self._name, st,en)
				cg = dg.g(id=subseq, font_size=size, display='none')
				group.add(cg)
				self.text(text="{:d}".format(st), x=x1, y=y1-self._domain_height*0.75, 
					group=cg, size=size, align='right')
				self.text(text="{:d}".format(en), x=x2, y=y1-self._domain_height*0.75,
					group=cg, size=size, align='left')

	def subsequences(self, domain_boundaries, group=None, size=9):
		dg = self._drawing
		if group is None:
			group = dg
		if len(domain_boundaries) > 0:
			for (st, en) in domain_boundaries:
				(x1, y1, x2, y2) = self.lineBounds(st, en)
				subseq = "{:s}-{:d}-{:d}-seq".format(self._name, st,en)
				cg = dg.g(id=subseq, display='none')
				group.add(cg)
				self.text(text="{:s}".format(self._seq[(st-1):en]), x=(x1+x2)/2.0, y=y1+self._domain_height, 
					group=cg, font='Lucida Console', size=size, align='center')




if __name__=='__main__':
	sd = SequenceDrawing(width=20, height=10)
	sd.cartoon(1, "AAAAAAAAAAGGGGGGGGGG", [(1,9),(11,19)], ['red','blue'])
	sd.sequence(2, "AAAAAAAAAAGGGGGGGGGG", [(1,9),(11,19)], ['red','blue'])
	sd.label(1, "Npl3")
	sd.write('test.svg')
