package com.zx.MultiThreaded;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

public class FileReader {

	private FileChannel fileChanne;
	private ByteBuffer byteBuffer;
	private int bufferSize;
	private long start;

	// File offset of byte[0] of the current byteBuffer content.
	// Updated every time a new buffer is read from the channel.
	private long bufferBasePos;

	// File offset of the first byte of the most recently returned line.
	// Callers use this instead of fileChannel.position() to detect overlap boundaries,
	// because fileChannel.position() reflects the pre-fetched channel position which
	// is always bufferSize bytes ahead of where the current line actually starts.
	private long lastLineStart;

	public FileReader(FileChannel fileChannel, int bufferSize, long start) throws IOException {
		this.fileChanne = fileChannel;
		this.bufferSize = bufferSize;
		this.start = start;
		fileChannel.position(start);
		this.bufferBasePos = start;
		this.lastLineStart = start;
	}

	/** Returns the file offset of the first byte of the most recently returned line. */
	public long getLastLineStart() {
		return lastLineStart;
	}

	public String readline() throws IOException {

		if (byteBuffer == null) {
			byteBuffer = ByteBuffer.allocate(bufferSize);

			int len = fileChanne.read(byteBuffer);

			if (len == -1) {
				return null;
			}

			byteBuffer.flip();
			// bufferBasePos is already set to start in the constructor
		}

		byte[] bb = new byte[bufferSize];
		int i = 0;

		// Capture the file offset of the first byte of THIS line.
		long thisLineStart = bufferBasePos + byteBuffer.position();

		while (true) {

			while (byteBuffer.hasRemaining()) {

				byte b = byteBuffer.get();

				if ('\n' == b || '\r' == b) {

					if (byteBuffer.hasRemaining()) {
						byte n = byteBuffer.get();

						if ('\n' != n) {
							byteBuffer.position(byteBuffer.position() - 1);
						}

					} else {

						byteBuffer.clear();
						// Record the base position of the next buffer before reading it.
						bufferBasePos = fileChanne.position();
						int len = fileChanne.read(byteBuffer);

						byteBuffer.flip();

						if (len != -1) {
							byte n = byteBuffer.get();

							if ('\n' != n) {
								byteBuffer.position(byteBuffer.position() - 1);
							}
						}

					}

					lastLineStart = thisLineStart;
					return new String(bb, 0, i);

				} else {

					if (i >= bb.length) {

						bb = Arrays.copyOf(bb, bb.length + bufferSize + 1);
					}

					bb[i++] = b;
				}

			}

			// Buffer exhausted — load the next chunk.
			byteBuffer.clear();
			bufferBasePos = fileChanne.position();
			int len = fileChanne.read(byteBuffer);
			byteBuffer.flip();

			if (len == -1 && i == 0) {
				return null;
			}

		}

	}

	public void close() throws IOException {
		this.fileChanne.close();
	}


}
