# Image Editor
#### Jun-Yu (Andrew) Chen
## Introduction
An Image Editor written in C++ implementing various effects on images. The editing program allows loading of one or more images and performs various operations on them. Similar to creating a miniature photoshop.
<p align="center">
  <img src="https://user-images.githubusercontent.com/64970325/151760829-2ea0928a-dc7c-4d5a-8369-e647ff898cf6.PNG" />
  <img src="https://user-images.githubusercontent.com/64970325/151760970-7a260139-7260-453f-9960-3a55a896b61c.PNG" />
</p>

## Building and Running
Can be built directly with CMake with the given CMakeLists text file, preferred building/running platform would be Visual Studio 2019 on Windows 10.

## Calling of Tasks Done
| Effect  | Input Keyword |
| ------------- | ------------- |
| Load Picture (In .tga format) | load \<source\> |
| Grayscale | gray |
| Uniform Quantization | quant-unif |
| Populosity | quant-pop |
| Naive Threshold Dithering | dither-thresh |
| Brightness Preserving Threshold Dithering | dither-bright |
| Random Dithering | dither-rand |
| Clustered Dithering | dither-cluster |
| Floyd-Steinberg Dithering | dither-fs |
| Box Filter | filter-box |
| Bartlett Filter | filter-bartlett |
| Gaussian Filter | filter-gauss |
| Arbitrary-Size Gaussian Filter | filter-gauss-n <size(odd int)> |
| Half Size | half |
| Double Size | double |
| Save | save \<filename\> |

## Implementation
Done as coursework for CS3026 Introduction to Computer Graphics in NTUST, Fall 2021. More details about the course project could be seen in [this link](http://dgmm.csie.ntust.edu.tw/?ac1=courprojdetail_CG2012F_3&id=5ecf7b7a5118c&sid=614a94d120553).

##### © 2022 Andrew (Jun-Yu) Chen
