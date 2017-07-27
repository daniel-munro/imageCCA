run("Bio-Formats Importer", "open=" + getArgument + " autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_3");
run("Set Measurements...", "  mean redirect=None decimal=0");
inc = round(1000 / 16);
c = 180;
goodX = newArray();
goodY = newArray();
x = 1;
y = 1;
while (x + 2 * inc < getWidth) {
	while (y + 2 * inc < getHeight) {
		makeRectangle(x, y, inc, inc);
		run("Measure");
		if (getResult('Mean', nResults - 1) < c) {
			makeRectangle(x, y - inc, inc, inc);
			run("Measure");
			u = getResult('Mean', nResults - 1);
			makeRectangle(x, y + inc, inc, inc);
			run("Measure");
			d = getResult('Mean', nResults - 1);
			makeRectangle(x - inc, y, inc, inc);
			run("Measure");
			l = getResult('Mean', nResults - 1);
			makeRectangle(x + inc, y, inc, inc);
			run("Measure");
			r = getResult('Mean', nResults - 1);
			if (u < c && d < c && l < c && r < c) {
				goodX = Array.concat(goodX, x);
				goodY = Array.concat(goodY, y);
			}
		}
		y += inc;
	}
	x += inc;
	y = 1;
}
run("Close");
print("good tiles: " + goodX.length);
if (goodX.length > 0) {
	choice = floor(random * goodX.length);
	x = round(goodX[choice] * 16);
	y = round(goodY[choice] * 16);

	run("Bio-Formats Importer", "open=" + getArgument + " autoscale color_mode=Default crop view=Hyperstack stack_order=XYCZT series_1 x_coordinate_1=" + x + " y_coordinate_1=" + y + " width_1=1000 height_1=1000");
	saveAs("Jpeg", "~/images/tmp.jpg");
}
