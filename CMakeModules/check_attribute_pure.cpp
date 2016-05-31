int pure_function() __attribute__ ((pure));

int pure_function() {
    return 0;
}

int main()
{
	return pure_function();
}
