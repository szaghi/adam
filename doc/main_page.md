project: ADAM
src_dir: ../src/app/nasto
exclude_dir: ../src/third_party/
             ../src/third_party_manual/
output_dir: html/publish/
project_github: https://github.com/szaghi/adam
summary: Adaptive Mesh Refinement (AMR) with Immersed Boundary (IB) fluid dynamic solver tailored for High Performance GPU Computing
author: Stefano Zaghi
github: https://github.com/szaghi
email: stefano.zaghi@gmail.com
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
               markdown.extensions.extra
docmark: <
display: public
         protected
         private
source: true
warn: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README-ADAM.md!}
