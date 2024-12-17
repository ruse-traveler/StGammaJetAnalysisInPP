TLegend* MakeTextToAddToClosureTestForPaper() {
  double xy[4] = {0.67, 0.92, 0.94, 0.965};
  TLegend* txt = new TLegend(xy[0], xy[1], xy[2], xy[3]);
  txt -> SetTextFont(42);
  txt -> SetTextSize(0.085);
  txt -> SetFillColor(0);
  txt -> SetLineColor(0);
  txt -> AddEnty((TObject*) 0, "p+p, #sqrt{s} = 200 GeV", "");
  return txt;
}