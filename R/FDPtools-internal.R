.cdev <-
function () 
{
    x11(type = "Xlib")
}
.xF <-
function () 
{
    source(system("xclip -o", intern = T))
}
.xR <-
function () 
{
    source("clipboard")
}
