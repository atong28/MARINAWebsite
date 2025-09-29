#!/bin/bash

# SPECTRE Docker Restart Script
# This script performs a hard restart of the Docker containers with logs

set -e

echo "🔄 SPECTRE Docker Restart Script"
echo "================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Docker is running
if ! docker info >/dev/null 2>&1; then
    print_error "Docker is not running. Please start Docker and try again."
    exit 1
fi

# Check if docker-compose is available
if ! command -v docker-compose >/dev/null 2>&1; then
    print_error "docker-compose is not installed. Please install docker-compose and try again."
    exit 1
fi

print_status "Stopping existing containers..."
docker-compose down --remove-orphans

print_status "Removing old images to force rebuild..."
docker-compose build --no-cache

print_status "Starting containers in detached mode..."
PORT=8001 docker-compose up -d

# Wait a moment for containers to start
sleep 3

# Check if container is running
if docker-compose ps | grep -q "Up"; then
    print_success "Containers started successfully!"
    
    # Show container status
    echo ""
    print_status "Container Status:"
    docker-compose ps
    
    echo ""
    print_status "Recent logs:"
    docker-compose logs --tail=20
    
    echo ""
    print_status "To follow logs in real-time, run:"
    echo -e "${YELLOW}docker-compose logs -f${NC}"
    
    echo ""
    print_status "To stop containers, run:"
    echo -e "${YELLOW}docker-compose down${NC}"
    
    # Ask if user wants to follow logs
    echo ""
    read -p "Would you like to follow logs now? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_status "Following logs (Ctrl+C to exit)..."
        docker-compose logs -f
    fi
    
else
    print_error "Failed to start containers. Checking logs..."
    docker-compose logs
    exit 1
fi
