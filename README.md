MCREEP
------

MCREEP package fits various pre-defined models to creep data.

The package can process both tensile and indentation creep measuremenents <br>
in such a way that the results from the two methods are comparable.

* The creep data:
    - A text file with two columns.
    - The 1st column = time, the 2nd column = deformation.
	- Deformation = strain (for tensile) or penetration depth (for indentation).
* The pre-defined models within the pagkage are:
    - Simple empirical models = Power law, Nutting's law.
    - Phenomenological models = elasto-visco-plastic models = EVP models.
* The theory and models are described in the open-acces publication:
	- *Materials*, submitted (minor revisions requested).
	- If you use MCREEP in your research, **please cite** the publication above.

Quick start
-----------
1. Install MCREEP: `pip install mcreep`
2. Download the [demo package](https://mirekslouf.github.io/mcreep/docs),
   unzip it, and follow the instructions in *readme* file.
3. Basically, you just run the unzipped PY-script and see the results.
	- The PY-script should run anywhere, but the recommended environment is
      [Spyder](https://www.spyder-ide.org).
	- The output (in Spyder) should look like in this 
      [PrintScreen](https://mirekslouf.github.io/mcreep/docs/images/mcreep_printscreen.png).

Installation, documentation and examples
--------------------------
* MCREEP is available at
  [PyPI](https://pypi.org/project/mcreep/)
  &rArr; installation: `pip install mcreep`
* Home page of the project is at
  [GitHub](https://github.com/mirekslouf/mcreep/) 
  and [GitHub pages](https://mirekslouf.github.io/mcreep/).
* [Docs and examples](https://mirekslouf.github.io/mcreep/docs/)
  are summarized at GitHub pages.
* All modules, classes and functions are documented by docstrings.
* [Full documentation](https://mirekslouf.github.io/mcreep/docs/pdoc.html/index.html)
  was auto-generated from docstrings by [PDoc](https://pdoc.dev).

Versions of MCREEP
------------------
* Version 1.0 = finalized, working for both tensile and indentation creep
* Version 1.1 = a few important improvements, better documentation
