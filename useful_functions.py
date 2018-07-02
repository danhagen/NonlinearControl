import numpy as np
import os.path
import matplotlib._pylab_helpers
from matplotlib.backends.backend_pdf import PdfPages

def return_length_of_nonzero_array(X):
	"""
	Takes in a numpy.ndarray X of shape (m,n) and returns the length of the array that removes any trailing zeros.
	"""
	assert str(type(X))=="<class 'numpy.ndarray'>", "X should be a numpy array"
	assert np.shape(X)[1]!=1, "X should be a wide rectangular array. (m,1) is a column, therefore a nonzero X of this shape will return 1 (trivial solution). Transpose X to properly identify nonzero array length."
	assert np.shape(X)!=(1,1), "Check input. Should not be of shape (1,1) (trivial solution)."
	if (X[:,1:]!=np.zeros(np.shape(X[:,1:]))).all():
		return(np.shape(X)[1])
	else:
		return(np.argmax((X[:,1:] == np.zeros(np.shape(X[:,1:]))).sum(axis=0) == np.shape(X[:,1:])[0])+1)


def save_figures(BaseFileName,**kwargs):
	figs = kwargs.get("figs",
		[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
		)

	i = 1
	FileName = BaseFileName + "_" + "{:0>2d}".format(i) + ".pdf"
	if os.path.exists(FileName) == True:
		while os.path.exists(FileName) == True:
			i += 1
			FileName = BaseFileName + "_" + "{:0>2d}".format(i) + ".pdf"
	PDFFile = PdfPages(FileName)
	if len(figs)==1:
		PDFFile.savefig(figs[0])
	else:
		[PDFFile.savefig(fig) for fig in figs]
	PDFFile.close()
