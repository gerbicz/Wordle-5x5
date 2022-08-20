# Wordle-5x5
Fast search for Matt Parker's 5 five letter words that has 25 different letters from a-z.
Single, and multi threaded codes.

The base times: 0.05 seconds on single thread, and 0.02 on 8 threads
for the original 12972 words list, see vocabulary.txt. There are 11 solutions.

On the larger vocabulary it needs roughly 0.1 sec to find all the 831 solutions on a single thread,
see words_alpha.txt at https://github.com/dwyl/english-words for this word list
(it will consider only the 5 letter words).
