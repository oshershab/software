# Use a base image with necessary build tools
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && \
    apt-get install -y build-essential gfortran wget python3 python3-pip

# Download and extract the file (assuming a different format or URL)
RUN wget https://www.nupack.org/download/nupack-4.0.tar && \
    tar -xf nupack-4.0.tar && \
    cd nupack-4.0 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf nupack-4.0 nupack-4.0.tar

# Set the working directory
WORKDIR /tool

# Copy server code and requirements.txt
COPY tool/server.py .
COPY requirements.txt .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Set the command to run your server
CMD ["python3", "server.py"]
