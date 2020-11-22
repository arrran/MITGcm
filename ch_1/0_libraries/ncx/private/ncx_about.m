function ncx_about

bg = get(gca,'color');
Message={
  'NCX  -  NetCdf eXplorer',
  'NetCDF Visualisation utility'
  '',
  'version 1.0 December 2005',
  '',
  'Created by:'
  'Martinho Marta Almeida, martinho@fis.ua.pt',
  'Physics Department',
  'Aveiro University'
  'Portugal',
  '',
  'http://neptuno.fis.ua.pt/~mma/ncx'
};
Title='NCX, about';
mb=msgbox(Message,Title,'help','modal');
set(mb,'color',bg);
