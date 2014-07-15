#pragma once
/*************/
// Region fill
namespace Regionfill
{

#define MAX 10000		/* max depth of stack */
#define pixelread(X,Y) ((Pixel)img[X][Y])
#define pixelwrite(X,Y,V) img[X][Y] = V
  typedef struct {		/* window: a discrete 2-D rectangle */
    int xmin, ymin;			/* xmin and ymin */
    int xmax, ymax;			/* xmax and ymax (inclusive) */
  } Window;

  typedef short int Pixel;
  typedef struct {short y, xl, xr, dy;} Segment;
  /*
   * Filled horizontal segment of scanline y for xl<=x<=xr.
   * Parent segment was on line y-dy.  dy=1 or -1
   */

  /*
   *  push a segment with row Y, min-column XL, max-column XR, and direction DY
   *  into a stack
   */
#define PUSH(Y, XL, XR, DY)	/* push new segment on stack */         \
      if (sp<stack+MAX && Y+(DY)>=win->ymin && Y+(DY)<=win->ymax)       \
      {sp->y = Y; sp->xl = XL; sp->xr = XR; sp->dy = DY; sp++;}

  /*
   *  pop a segment with row Y+DY, min-column XL, max-column XR
   */
#define POP(Y, XL, XR, DY)	/* pop segment off stack */     \
  {sp--; Y = sp->y+(DY = sp->dy); XL = sp->xl; XR = sp->xr;}

  void regionfill(unsigned char **img, int xsize, int ysize, int x, int y, unsigned char nnv)
  {
    int l, x1, x2, dy;
    Pixel  ov, nv;	/* old pixel value */
    Segment stack[MAX], *sp = stack;	/* stack of filled segments */
    Window *win, wwin;	/* screen window */
	
    nv = (Pixel)nnv;
    wwin.xmin = 0;  wwin.xmax = xsize-1;
    wwin.ymin = 0;  wwin.ymax = ysize-1;
    win = &wwin;
    if (x<win->xmin || x>win->xmax || y<win->ymin || y>win->ymax || (ov=pixelread(x,y))==nv) return;
    PUSH(y, x, x, 1);			/* needed in some cases */
    PUSH(y+1, x, x, -1);		/* seed segment (popped 1st) */
	
    while (sp>stack) {
      /* pop segment off stack and fill a neighboring scan line */
      POP(y, x1, x2, dy);
      /*
       * segment of scan line y-dy for x1<=x<=x2 was previously filled,
       * now explore adjacent pixels in scan line y
       */
      for (x=x1; x>=win->xmin && pixelread(x, y)==ov; x--)
      {
        pixelwrite(x, y, nv);
      }

			
      if (x>=x1) goto skip;
      l = x+1;
      if (l<x1) PUSH(y, l, x1-1, -dy);		/* leak on left? */
      x = x1+1;
      do {
        for (; x<=win->xmax && pixelread(x, y)==ov; x++)
        {
          pixelwrite(x, y, nv);
        }
				
        PUSH(y, l, x-1, dy);
        if (x>x2+1) PUSH(y, x2+1, x-1, -dy);	/* leak on right? */
     skip:	    for (x++; x<=x2 && pixelread(x, y)!=ov; x++)
        {
        };
        l = x;
      } while (x<=x2);
    }
  }
  /***********************/
}
