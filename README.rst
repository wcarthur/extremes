extremes - a package for fitting extreme value distributions
============================================================

This package provides users with ways to interrogate the relational
database created by the Tropical Cyclone Risk Model (TCRM -
http://github.com/GeoscienceAustralia/tcrm) to fit return period
curves using the Generalised Pareto Distribution (GPD). 

The codes here enable a user to extract all simulated wind speeds for
a location and use two different methods for fitting the GPD. The
first method is the `fit` method available in all
`scipy.stats.rv_continuous` distributions, which uses maximum
likelihood estimates for the fitting. The second method is described
by `Sanabria &
Cechet <https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&uact=8&ved=0ahUKEwidmJn0lYPWAhUKwLwKHX22Ah0QFggxMAM&url=https%3A%2F%2Fd28rz98at9flks.cloudfront.net%2F65052%2FRec2007_012.pdf&usg=AFQjCNFKAeLHHq0eKZ7Q--9zm6sxWXWJuA>`_
(2007), which applies an iterative sampling technique to determine the
optimum threshold parameter for the GPD.

Using `extremes`
----------------

Dependencies:
~~~~~~~~~~~~~

- TCRM (v2.0)
- matplotlib
- scipy
- seaborn



You will need to install the TCRM package to access the `Utilities` module, which contains functions for configuration, log files, parallel processing and versioning. 

Contributing:
-------------

Anyone is welcome to contribute to the package. 

1. Fork the repository - you will need to have a GitHub account to do
   this. Click the 'Fork' button near the top right of the page. This
   will create a copy of the code under your account on the GitHub
   server.

2. Clone the repository to your local disk (using SSH authentication)::

     $ git clone git@github.com:YourLogin/extremes.git
     $ cd extremes

3. Create a new branch to hold your changes and start making changes::
     
     $ git checkout -b my-feature

4. Commit your changes to your fork of the repository::
     
     $ git add <modified_files>
     $ git commit 

   Then push them to GitHub with ::
     
     $ git push -u origin my-feature

5. Go to the webpage of your fork and click "Pull Request" to send
   your changes to me.

Please try to follow these guidelines when contributing new code. 

- All public methods should have informative docstrings with parameter
  definitions, sample usage, etc.

- Where possible, include a paragraph of narrative documentation,
  including links to references in the literature where possible (PDF
  links!)

- New functions should include unit tests that can be picked up using nose

- Check the code with pyflakes, pylint or some other PEP8 checker
  (autopep8 is a quick way to fix redundant PEP8 errors)


For more details on using Git, check out the `Git documentation
<http://git-scm.com/documentation>`_.

Contact: Craig Arthur, craig.arthur@ga.gov.au

