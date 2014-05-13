import numpy, pyfits, types, itertools
def mwrfits(filename, arraylist, namelist=None, header=None):
	""" 
	Writes the list of numpy.arrays arraylist as a FITS table filename
	using namelist as list of names. 
	Arraylist can be dictionary with arrays as values and names as keys. 
	Also Arraylist can be numpy-record-array.

	Example:
	mwrfits('/tmp/xx.fits',[arr,arr1],['X','Y'])
	Or :
	mwrfits('test.fits',{'X':arr,'Y':arr1})
	Or:
	data = numpy.zeros((4,),dtype=[('run','i4'),('rerun','f8'),('zz','b')])
	mwfits('test1.fits',data)
	
	Keep in mind that when you used a dictionary, the order of columns in the
	fits file is not guaranteed
	"""

	tmplist=[]
	print(type(arraylist))
	print('start')
	if isinstance(arraylist,numpy.ndarray):
		print('1')

		if arraylist.dtype.type is numpy.void:
			print('1.2')
			iter=itertools.izip(arraylist.dtype.names, itertools.imap (arraylist.__getitem__ , arraylist.dtype.names))
		else:
		#NOTE TO SELF: This else statement has been added in... not sure what type of iter is needed though
			iter = numpy.ndenumerate(arraylist)
			print('1.3')
			
	else:
		print('2')
		if isinstance(arraylist,types.ListType):
			print('2.2')
			iter= zip(namelist, arraylist)
		elif isinstance(arraylist,types.DictType):
			print('2.3')
			iter= arraylist.iteritems()

	for name, arr in iter:
		if arr.dtype.type==numpy.int8:
			format='I'
		elif arr.dtype.type==numpy.int16:
			format='I'
		elif arr.dtype.type==numpy.int32:
			format='J'
		elif arr.dtype.type==numpy.int64:
			format='K'
		elif arr.dtype.type==numpy.float32:
			format='E'
		elif arr.dtype.type==numpy.float64:
			format='D'
		elif arr.dtype.type==numpy.string_:
			format='%dA'%arr.dtype.itemsize
		else:
			raise Exception("Oops unknown datatype %s"%arr.dtype)
		test = numpy.array([arr])
		tmplist.append(pyfits.Column(name=name, array=test, format=format))
		
		#tmplist.append(pyfits.Column(name=name, format=format))

	#NOTE TO SELF: This line is the one throwing the 'tuple index out of range' error. Error caused by table.py dim = arr.shape[0]
	hdu = pyfits.new_table(tmplist)
	hdu.writeto(filename,clobber=True)