#ifndef ARENA_H
#define ARENA_H

#include <algorithm>
#include <iostream>
#include <list>

#define L1_CACHE_LINE_SIZE 64

static void* alloc_aligned(size_t size) {
#if defined(PLATFORM_MACOS)
	void *ptr;
    if (posix_memalign(&ptr, L1_CACHE_LINE_SIZE, size) != 0)
        ptr = nullptr;
    return ptr;
#elif defined(PLATFORM_WINDOWS)
	return _aligned_malloc(size, L1_CACHE_LINE_SIZE);
#endif
}

template <typename T>
static T* alloc_aligned(size_t count) {
	return (T*) alloc_aligned(count * sizeof(T));
}

static void free_aligned(void * ptr) {
#if defined(PLATFORM_MACOS)
	free(ptr);
#elif defined(PLATFORM_WINDOWS)
	_aligned_free(ptr);
#endif
}

class mem_arena {
public:
	mem_arena(size_t block_size = 262144) : block_size(block_size) {};
	~mem_arena() {
		// std::cerr << "destroying arena with " << total_allocated() << " bytes allocated.\n";
		free_aligned(current_block);
		for (auto& block : used_blocks)
			free_aligned(block.second);
		for (auto& block : available_blocks)
			free_aligned(block.second);
	}

	void* alloc(size_t n_bytes) {		
		n_bytes = (n_bytes + 15) & (~15);

		if (current_pos + n_bytes > current_alloc_size) {
			if (current_block) {
				used_blocks.emplace_back(current_alloc_size, current_block);
				current_block = nullptr;
			}

			for (auto iter = available_blocks.begin(); iter != available_blocks.end(); ++iter) {
				if (iter->first >= n_bytes) {
					current_alloc_size = iter->first;
					current_block = iter->second;
					available_blocks.erase(iter);
					break;
				}
			}

			if (!current_block) {
				current_alloc_size = std::max(n_bytes, block_size);
				current_block = alloc_aligned<unsigned char>(current_alloc_size);
			}
			current_pos = 0;
		}
		void *ret = current_block + current_pos;
		current_pos += n_bytes;
		return ret;
	};

	template <typename T, typename... Args>
	T* alloc(Args&&... args) {
		T* ret = (T*) alloc(sizeof(T));
		new (ret) T(std::forward<Args>(args)...);
		return ret;
	};
	
	void reset() {
		current_pos = 0;
		available_blocks.splice(available_blocks.begin(), used_blocks);
	}

	size_t total_allocated() const {
		size_t total = current_alloc_size;
		for (const auto &block : used_blocks)
			total += block.first;
		for (const auto &block : available_blocks)
			total += block.first;
		return total;
	}

private:
	size_t block_size;
	size_t current_pos = 0, current_alloc_size = 0;
	unsigned char *current_block = nullptr;
	std::list<std::pair<size_t, unsigned char*>> used_blocks, available_blocks;
};

#endif