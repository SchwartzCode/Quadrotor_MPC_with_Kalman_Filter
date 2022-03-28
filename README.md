# MPC Control of a Quadrotor with Wind Disturbance

## README

### README

#### README

##### README

readme

Things to do:
- edit functions in dynamis.jl to match the dynamics laid out in main_script.ipynb (Corinne)
- find reasonable values for quadrotor parameters (both)
- make a quadrotor struct to store the above params? (both?)
- implement toggle-able wind disturbance in dynamics (Corinne)\
    - [brown note](https://www.youtube.com/watch?v=mQFL-NLh0O8)
    - [okay but actually this sounds really nice](https://www.youtube.com/watch?v=hXetO_bYcMo)
- some sort of test to show dynamics model works - maybe simple LQR controller? (both)
- write out constraints for MPC (Jonathan)
- re-write MPC functions to work for this system (Jonathan)

Some notes:
- to add a package, you need to add it to Project.toml (and probably delete Manifest.toml and rebuild)
