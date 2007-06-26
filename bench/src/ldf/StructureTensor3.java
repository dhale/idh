/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ldf;

/*
StructureTensor3
  CompressedStructureTensor3
    construct from six 3-D arrays g1g1, g1g2, g1g3, g2g2, ...
    use Denning's compression for normal vectors
    quantize eigenvectors u and w; get v from cross-product
    u 16 bits 
    w 16 bits 
    eu 32 bits
    ev/eu 8 bits
    ew/eu 8 bits
    total is 80 bits = 2.5 floats/sample
    use coarser grid to store filters
      store barycentric weights on finer grid
      need to be able to store information at all levels
      
    float[] getU(int i1, int i2, int i3)
    float[] getE(int 
  BlockedStructureTensor3
    average eight nearest structure tensors

SphereSampling

*/
/**
 * Interface for 3-D structure tensors.
 * @author Dave Hale, Colorado School of Mines
 * @version 2007.06.26
 */
public interface StructureTensor3 {
}
