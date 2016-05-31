int const_function() __attribute__ ((const));

int const_function() {
    return 0;
}

int main()
{
	return const_function();
}
