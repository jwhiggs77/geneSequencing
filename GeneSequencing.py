#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

	def __init__(self):
		pass

	# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
	# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
	# how many base pairs to use in computing the alignment

	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		# w1 is seqence at the top of the table
		w1 = " " + seq1[:align_length]
		# w2 is seqence at the side of the table
		w2 = " " + seq2[:align_length]

		if seq1[:align_length] == seq2[:align_length]:
			score = MATCH * (len(w1) - 1)
			alignment1 = seq1[:100]
			alignment2 = seq2[:100]
			return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
		# IF NOT BANNED
		if banded == False:
			# make table
			rows, cols = (len(w1), len(w2))
			newCell = (0, 0)
			table = [[newCell for i in range(cols)] for j in range(rows)]

			# initialize the table
			for i in range(1, rows):
				cost = INDEL * i
				backPnt = {"ri": i - 1, "ci": 0}
				cell = (cost, backPnt)
				table[i][0] = cell
			for i in range(1, cols):
				cost = INDEL * i
				backPnt = {"ri": 0, "ci": i - 1}
				cell = (cost, backPnt)
				table[0][i] = cell

			for i in range(1, rows):
				for j in range(1, cols):
					# cost of left
					lowerCost = table[i][j - 1][0] + INDEL
					backPnt = {"ri": i, "ci": j - 1}
					# cost of above
					cost = table[i - 1][j][0] + INDEL
					if lowerCost > cost:
						lowerCost = cost
						backPnt = {"ri": i - 1, "ci": j}
					if w1[i] == w2[j]:
						# cost of match
						cost = table[i - 1][j - 1][0] + MATCH
						if lowerCost > cost:
							lowerCost = cost
							backPnt = {"ri": i - 1, "ci": j - 1}
					else:
						# cost of corner
						cost = table[i - 1][j - 1][0] + SUB
						if lowerCost > cost:
							lowerCost = cost
							backPnt = {"ri": i - 1, "ci": j - 1}
					# add cell to table
					cell = (lowerCost, backPnt)
					table[i][j] = cell
			alignment1, alignment2 = self.makeUnbannedAlign(table, w1, w2)
			score = table[rows - 1][cols - 1][0]
		# IF BANNED
		else:
			if abs(len(seq1) - len(seq2)) > 500:
				score = math.inf
				alignment1 = "No Alignment Possible"
				alignment2 = "No Alignment Possible"
				return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
			# make table
			rows, cols = (len(w1), len(w2))
			if rows > cols:
				largerWordSize = rows
			else:
				largerWordSize = cols
			table = []
			for i in range(rows):
				if i < 3:
					d = i + 1
				elif (largerWordSize - i) < 4:
					d -= 1
				else:
					d = 4
				band = 3 + d
				col = []
				for j in range(band):
					col.append(0)
				table.append(col)

			#initialize
			table[0][0] = newCell = (0, 0)
			band = 4
			for i in range(1, band):
				cost = INDEL * i
				backPnt = {"ri": i - 1, "ci": 0}
				cell = (cost, backPnt)
				table[i][0] = cell
			for i in range(1, band):
				cost = INDEL * i
				backPnt = {"ri": 0, "ci": i - 1}
				cell = (cost, backPnt)
				table[0][i] = cell

			shiftJ = 0
			for i in range(1, rows):
				band = len(table[i])
				if (i > 3):
					shiftJ += 1
				for j in range(band):
					if (i > 3):
						jShf = j + 1 # j shifted
						jm = j + shiftJ
					else:
						jShf = j
						jm = j
					# cost of left
					if j > 0:
						cost = table[i][j - 1][0] + INDEL
						lowerCost = cost
						backPnt = {"ri": i, "ci": j - 1}
					else:
						lowerCost = math.inf
					# cost of above
					if jShf < len(table[i - 1]):
						cost = table[i - 1][jShf][0] + INDEL
						if lowerCost > cost:
							lowerCost = cost
							backPnt = {"ri": i - 1, "ci": jShf}
					# if there is a match
					if w1[i] == w2[jm]:
						# cost of match
						cost = table[i - 1][jShf - 1][0] + MATCH
						if lowerCost > cost:
							lowerCost = cost
							backPnt = {"ri": i - 1, "ci": jShf - 1}
					else:
						# cost of corner
						cost = table[i - 1][jShf - 1][0] + SUB
						if lowerCost > cost:
							lowerCost = cost
							backPnt = {"ri": i - 1, "ci": jShf - 1}
					# add cell to table
					cell = (lowerCost, backPnt)
					table[i][j] = cell
			lastCol = len(table[rows - 1]) - 1
			score = table[rows - 1][lastCol][0]
			alignment1, alignment2 = self.makeBannedAlign(table, w1, w2)

		return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}


	def makeUnbannedAlign(self, table, w1, w2):
		align1 = ''
		align2 = ''
		# create alignments
		END = {"ri": 0, "ci": 0}
		i = len(w1) - 1
		j = len(w2) - 1
		curCell = {"ri": i, "ci": j}
		pointer = table[i][j][1]
		while curCell != END:
			rIndex = pointer.get("ri")
			cIndex = pointer.get("ci")
			if rIndex < curCell.get("ri") and cIndex < curCell.get("ci"):  # if points diagonal
				align1 += w1[i]
				align2 += w2[j]
				i -= 1
				j -= 1
			elif cIndex < curCell.get("ci"):  # if points left
				align1 += "-"
				align2 += w2[j]
				j -= 1
			elif rIndex < curCell.get("ri"):  # if points above
				align1 += w1[i]
				align2 += "-"
				i -= 1
			curCell = {"ri": rIndex, "ci": cIndex}
			pointer = table[rIndex][cIndex][1]
		if len(align1) > 100 or len(align2) > 100:
			align1 = align1[-100:]
			align2 = align2[-100:]
		align1 = ''.join(reversed(align1))
		align2 = ''.join(reversed(align2))
		return align1, align2

	def makeBannedAlign(self, table, w1, w2):
		align1 = ''
		align2 = ''
		# create alignments
		END = {"ri": 0, "ci": 0}
		i = len(w1) - 1
		j = len(w2) - 1
		curCell = {"ri": i, "ci": len(table[i]) - 1}
		pointer = table[i][len(table[i]) - 1][1]
		while curCell != END:
			rIndex = pointer.get("ri")
			cIndex = pointer.get("ci")
			if rIndex < curCell.get("ri") and cIndex > curCell.get("ci"):  # if points above
				align1 += w1[i]
				align2 += "-"
				i -= 1
			elif cIndex < curCell.get("ci") and rIndex >= curCell.get("ri"):  # if points left
				align1 += "-"
				align2 += w2[j]
				j -= 1
			elif rIndex < curCell.get("ri"):  # if points diagonal
				align1 += w1[i]
				align2 += w2[j]
				i -= 1
				j -= 1
			curCell = {"ri": rIndex, "ci": cIndex}
			pointer = table[rIndex][cIndex][1]
		if len(align1) > 100 or len(align2) > 100:
			align1 = align1[-100:]
			align2 = align2[-100:]
		align1 = ''.join(reversed(align1))
		align2 = ''.join(reversed(align2))
		return align1, align2