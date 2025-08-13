def convergence():
    return 1e-5

def scfmaxcycles():
    return 1000

def stylesheetqt(modstyle):
    styleProgressBar = """
QWidget:default
{
    font:bold;
}
QProgressBar
{
    
    border-radius: 5px;
    text-align: center;
}
"""
    styleData = """
QWidget
{
    color: white;
    font: bold;
}
QProgressBar
{
    border: 1px solid black;
    border-radius: 5px;
    text-align: center;
}
QProgressBar::chunk
{
    background-color: black;
}
"""
    if modstyle == "styleProgressBar":
        return styleProgressBar
    elif modstyle == "styleData":
        return styleData
