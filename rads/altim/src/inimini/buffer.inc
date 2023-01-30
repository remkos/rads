*     Defines a buffer for file operations

* nbytes	number of bytes for the buffer
* bufsiz	number of observations for the buffer
* minobs	lower observation number in current buffer
* maxobs	higher observation number in current buffer
* fd		file descriptor
* block		current block of nbytes loaded in buffer

      integer nbytes,bufsiz,minobs,maxobs,fd,block
      parameter (nbytes=320000)
      integer*1 buffer(nbytes)
      common /buffers/ buffer,bufsiz,minobs,maxobs,fd,block
