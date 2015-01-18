Introduction to AssayTools
--------------------------

Start by loading up an example assay dataset from disk. AssayTools will automatically parse the file extension and use the appropriate loader. ::

  >>> import assaytools as at
  >>> assay = at.load('assay.xml')
  >>> print assay

