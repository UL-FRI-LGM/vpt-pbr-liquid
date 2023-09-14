
# Physically based rendering of liquid environments

A repository containing code developed as part of the Bachelor thesis at the University of Ljubljana, Faculty of Computer and Information Science by [Aljaž Šmajcelj](https://github.com/asmaljcelj).

Related repository:
[https://github.com/asmaljcelj/volume_generation_and_simulation]

## Abstract
Rendering of volumetric data plays a vital role in modern society, yet it is far from an easy task. Frequently we desire physically correct depictions in which light refraction is also included when we are dealing with fluid inside a volume. In this work we present a way to generate and simulate a fluid inside a volume to clear up the fluid. The volume consists of multiple materials with different densities and other physical properties. We use this newly created volume in an already existing volumetric rendering framework with an extension that takes light refraction into account. For clear visualization of volumetric data, we created an automatic generation of transfer function to visually emphasize parts of volume with similar physical properties. We present the results in an existing web application of the framework where we evaluate these implemented extensions.


The thesis is accessible at:
[https://repozitorij.uni-lj.si/Dokument.php?id=147577]

Bibtex:
```
@thesis{Smajcelj2021,
    title = {Physically based rendering of liquid environments},
    author = {Šmajcelj, Aljaž},
    year = {2021},
    school = {University of Ljubljana, Faculty of Computer and Information Science},
    type = {Bachelor thesis},
    address = {Ljubljana, Slovenia},
    note = {Mentor: Ciril Bohak, Language: Slovenian, Slovenian title: Fizikalno osnovano upodabljanje tekočinskih okolij},
    url = {https://repozitorij.uni-lj.si/Dokument.php?id=147577}
}
```


...

Built on top of:

# [VPT: The Volumetric Path Tracing Framework](http://lgm.fri.uni-lj.si/portfolio-view/volumetric-path-tracing-framework/)

VPT is a volumetric path tracing framework targeted towards interactive
real-time data exploration. It works in both desktop and mobile environments.
It is built on top of WebGL 2 with no external dependencies.

![VPT](src/images/screenshot.jpg)

Visit the [portfolio page](http://lgm.fri.uni-lj.si/portfolio-view/volumetric-path-tracing-framework/) for more information.

## Building and running

You need only `node` to build the framework and to run it.

```bash
bin/packer
bin/server-node
```

There's a working build with a few demo datasets available [here](http://lgm.fri.uni-lj.si/~ziga).

## License

This project is licensed under the **GNU General Public License, version 3**.
See [LICENSE](LICENSE) for details.

