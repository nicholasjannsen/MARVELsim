Readme file
===========

Setting up a SPHINX documentation
---------------------------------

The following items is a guide:

* See the [documentation](https://www.sphinx-doc.org/en/master/index.html)
* Install nice theme [here](https://sphinx-themes.org/sample-sites/sphinx-rtd-theme/).

  
Troubleshooting
---------------

.. admonition:: rtd-theme      

   Note that sphinx 3.1 is required for the [markdown integration](https://www.sphinx-doc.org/en/master/usage/markdown.html) to work with the theme ``sphinx-rtd-theme``. However, note that there is a bug with sphinx-rtd-theme for the required sphinx version to use markdown (hence we here use ``rst`` format).


.. admonition:: rtd-theme
		
   To keep the nice layout of the ``rtd-theme`` you need to downgrade sphinx to version 1.8:

   .. code-block:: shell
   
      poetry add sphinx="1.8"
      

.. admonition:: docutils

   The package docutils=0.17 has a flaw for making itemized and enumarated lists in Sphinx. Downgrade to docutils=0.16 using:

   .. code-block:: shell

      poetry add docutils="0.16"
