from casatasks import flagdata
from pathlib import Path
import click
 
@click.command()
@click.argument("ms", type=Path)
def main(ms):
	flagdata(
		vis=ms, 
		mode="tfcrop",
		)
	pass

if __name__ == "__main__":
	main()
